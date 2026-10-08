"""Traceable comparison tables and monochrome figures from completed runs only."""

import csv
import json
import os
from pathlib import Path
import numpy as np
from .cli import write_json


def load_times(path):
    with path.open(newline="") as handle:
        return np.fromiter((float(r["dt_s"]) for r in csv.DictReader(handle)), dtype=float)


def compare_results(root, plot=True):
    root = Path(root)
    manifest = json.loads((root / "manifest.json").read_text())
    summaries = {}
    excluded = []
    for entry in manifest["scenarios"]:
        path = root / entry["name"] / "summary.json"
        summary = json.loads(path.read_text()) if path.exists() else {}
        if entry["status"] != "complete" or summary.get("status") != "complete":
            excluded.append(entry["name"])
            continue
        summaries[entry["name"]] = summary
    destination = root / "comparacao"
    destination.mkdir(exist_ok=True)
    rows = []
    for name, s in summaries.items():
        row = {"scenario": name, "sigma_s_eV": s["parameters"]["sigma_s"], "sigma_w_eV": s["parameters"]["sigma_w"],
               "nominal_sigma_E_eV": s["parameters"]["untruncated_sigma_E_eV"],
               "total_time_s": s["total_time_s"], "duration_s": s["duration_s"],
               "unique_wells": s["landscape"]["unique_wells"], "unique_saddles": s["landscape"]["unique_saddles"],
               "landscape_saddle_draws": s["landscape"]["saddle_draws"],
               "landscape_saddle_rejections": s["landscape"]["saddle_rejections"],
               "landscape_rejection_rate": s["landscape"]["rejection_rate"],
               "waiting_mean_s": s["waiting_times"]["mean_s"], "waiting_median_s": s["waiting_times"]["median_s"],
               "waiting_std_s": s["waiting_times"]["std_s"], "waiting_zero_count": s["waiting_times"]["zero_count"],
               "ones_fraction": s["digitization"]["ones_fraction"], "ties": s["digitization"]["ties"],
               "runs_count": s["nist_sp800_22_subset"]["runs_count"],
               "block_size": s["nist_sp800_22_subset"]["block_size"],
               "complete_blocks": s["nist_sp800_22_subset"]["complete_blocks"],
               "residual_bits": s["nist_sp800_22_subset"]["residual_bits"],
               "runs_prerequisite_met": s["nist_sp800_22_subset"]["runs_prerequisite_met"]}
        for sample, key in (("preflight", "preflight"), ("offered", "realized_offered_barriers"), ("selected", "selected_barriers")):
            for field in ("count", "mean_eV", "std_eV", "min_eV", "max_eV"):
                row[f"{sample}_{field}"] = s[key][field]
        for field in ("saddle_draws", "saddle_rejections", "rejection_rate"):
            row[f"preflight_{field}"] = s["preflight"][field]
        row.update({f"waiting_q{q.replace('%', '')}_s": v for q, v in s["waiting_times"]["quantiles_s"].items()})
        for test, result in s["nist_sp800_22_subset"]["tests"].items():
            row[f"{test}_statistic"] = result["statistic"]
            row[f"{test}_p"] = result["p_value"]
            row[f"{test}_status"] = result["status"]
        rows.append(row)
    with (destination / "metricas.csv").open("w", newline="") as handle:
        if rows:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    comparison = {"status": "complete" if not excluded else "incomplete", "excluded_scenarios": excluded,
                  "source_summaries": {name: f"../{name}/summary.json" for name in summaries}, "metrics": rows,
                  "limitations": ["one run per scenario; descriptive comparison", "median fitted on tested data",
                                  "three tests only, not the full NIST suite",
                                  "heterogeneity of saddles and wells varied jointly", "seeded numerical simulation"]}
    if plot and summaries:
        comparison["figures"] = comparison_figures(root, summaries, destination / "figuras")
    write_json(destination / "resumo.json", comparison)
    return comparison



def comparison_figures(root, summaries, output):
    """Four main comparisons, with the actual block size recorded in each summary."""
    output.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/apeiron-matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.facecolor": "white", "figure.facecolor": "white", "savefig.facecolor": "white"})
    markers = ["o", "s", "^", "D"]
    names = list(summaries)
    times = {name: load_times(root / name / "waiting_times.csv") for name in names}
    bits = {name: np.frombuffer((root / name / "bits.txt").read_bytes().strip(), dtype=np.uint8)-48 for name in names}
    records = []
    def save(fig, name, caption, **metadata):
        fig.savefig(output / f"{name}.png", dpi=350, bbox_inches="tight")
        fig.savefig(output / f"{name}.svg", bbox_inches="tight")
        plt.close(fig)
        records.append({"name": name, "caption": caption, **metadata})

    span = min(2000, min(len(v) for v in times.values()))
    fig, axes = plt.subplots(len(names), 1, figsize=(7, 2.25*len(names)), sharex=True, sharey=True, constrained_layout=True)
    for ax, name in zip(np.atleast_1d(axes), names):
        dt = times[name][:span]*1e12
        threshold = summaries[name]["digitization"]["median_dt_s"]*1e12
        ax.plot(np.arange(1,span+1),np.where(dt>0,dt,np.nan),color="black",linewidth=.5)
        if threshold > 0:
            ax.axhline(threshold,color="black",linestyle="--",label=f"Mediana {threshold:.3f} ps")
        ax.set(title=name,yscale="log",ylabel="Espera (ps)")
        ax.legend(fontsize=8)
    np.atleast_1d(axes)[-1].set_xlabel("Índice do salto")
    save(fig,"01-serie-temporal",f"Primeiros {span} saltos dos três cenários com escalas iguais e mediana de cada rodada. Zeros são omitidos somente do eixo logarítmico.",first_event=1,last_event=span)

    fractions = {}
    for name in names:
        tests=summaries[name]["nist_sp800_22_subset"]
        m=tests["block_size"]
        n=len(bits[name])//m
        fractions[name]=bits[name][:n*m].reshape(n,m).mean(axis=1) if n else np.array([])
    deviation=max([float(np.max(np.abs(v-.5))) for v in fractions.values() if v.size]+[.005])
    half_width=max(.01,deviation*1.2)
    fig,axes=plt.subplots(len(names),1,figsize=(7,2.1*len(names)),sharex=True,sharey=True,constrained_layout=True)
    for ax,name in zip(np.atleast_1d(axes),names):
        tests=summaries[name]["nist_sp800_22_subset"]
        values=fractions[name]
        ax.plot(np.arange(1,len(values)+1),values,color="black",marker="o",markersize=3,linewidth=.8)
        ax.axhline(.5,color="black",linestyle="--")
        ax.set(title=f"{name} | {tests['complete_blocks']} blocos de {tests['block_size']:,} bits",ylabel="Fração de uns",ylim=(max(0,.5-half_width),min(1,.5+half_width)))
    np.atleast_1d(axes)[-1].set_xlabel("Índice do bloco na sequência")
    save(fig,"02-equilibrio-blocos","Fração de uns nos blocos efetivamente usados no teste de frequência em blocos. Escalas verticais iguais e ampliadas em torno de 0,5.",block_sizes={name:summaries[name]["nist_sp800_22_subset"]["block_size"] for name in names})

    fig,ax=plt.subplots(figsize=(7,4.1),constrained_layout=True)
    counts=[summaries[name]["nist_sp800_22_subset"]["runs_count"] for name in names]
    reference=[2*len(bits[name])*bits[name].mean()*(1-bits[name].mean()) for name in names]
    upper=max([*counts,*reference]) or 1
    bars=ax.bar(names,counts,color="white",edgecolor="black",hatch=[["", "///", "..."][i % 3] for i in range(len(names))])
    for bar,count in zip(bars,counts):
        ax.text(bar.get_x()+bar.get_width()/2,count+.02*upper,f"{count:,}".replace(",","."),ha="center",fontsize=10)
    if len(set(reference))==1:
        ax.axhline(reference[0],color="black",linestyle="--",label="Referência do teste para o equilíbrio observado")
    else:
        ax.scatter(names,reference,color="black",marker="x",label="Referência do teste")
    ax.set(title="Trechos consecutivos de bits iguais",ylabel="Número de runs",ylim=(0,upper*1.16))
    ax.legend(fontsize=8,loc="lower left")
    save(fig,"03-runs","Contagens de trechos consecutivos de bits iguais e referência utilizada pelo teste de runs. Os totais iguais de zeros e uns não fixam a ordem dos bits.")

    fig,axes=plt.subplots(1,3,figsize=(7.5,4.1),sharey=True,constrained_layout=True)
    for ax,test,title in zip(axes,["monobit","block_frequency","runs"],["Frequência global","Frequência em blocos","Runs"]):
        for i,name in enumerate(names):
            result=summaries[name]["nist_sp800_22_subset"]["tests"][test]
            p=result["p_value"]
            if p is None:
                ax.text(i,.08,"Não aplicável",ha="center",rotation=90,fontsize=8)
            else:
                ax.plot(i,p,color="black",marker=markers[i%len(markers)],markerfacecolor="white",markersize=6)
                label="0 numérico" if p==0 else (f"{p:.2e}" if p<.001 else f"{p:.3f}")
                ax.annotate(label,(i,p),xytext=(0,9),textcoords="offset points",ha="center",fontsize=8,rotation=90 if 0<p<.001 else 0)
        ax.axhline(.01,color="black",linestyle="--",label="α = 0,01")
        ax.set(title=title,xticks=range(len(names)),xticklabels=names,xlim=(-.5,len(names)-.5),ylim=(-.04,1.14))
    axes[0].set_ylabel("Valor-p")
    axes[-1].legend(fontsize=8,loc="upper right")
    save(fig,"04-valores-p","Valores-p dos três testes com limite de rejeição de 0,01. A escala é linear e os valores pequenos são escritos junto aos marcadores. Zero numérico não é substituído por um número inventado.")
    write_json(output / "captions.json",records)
    return records
