"""Assemble manually authored chapters and actual tables; never generate conclusions."""
import argparse
import json
from pathlib import Path
import shutil

THESIS_ROOT = Path(__file__).resolve().parents[2] / 'tcc-apeiron-words'


def pt(value, digits=4):
    return f'{value:.{digits}f}'.replace('.',',')


def p_value(value):
    if value is None:
        return 'Não aplicável'
    if value == 0:
        return '0 numérico'
    if value >= .0001:
        return pt(value, 4)
    base, exponent = f'{value:.2e}'.split('e')
    superscript = str(int(exponent)).translate(str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))
    return base.replace('.', ',') + ' × 10' + superscript


def assemble(root, output, thesis_dir=None):
    original = (Path(thesis_dir) if thesis_dir is not None else THESIS_ROOT) / 'gerados'
    comparison=root/'comparacao'
    manifest=json.loads((root/'manifest.json').read_text())
    if manifest['status']!='complete': raise ValueError('cannot write thesis from incomplete runs')
    if (comparison/'reavaliacao-nist.json').exists():
        if json.loads((comparison/'reavaliacao-nist.json').read_text())['status'] != 'complete':
            raise ValueError('statistical reanalysis is incomplete')
    runs={e['name']:json.loads((root/e['name']/'summary.json').read_text()) for e in manifest['scenarios']}
    output.mkdir(parents=True,exist_ok=True)
    images=output/'imagens'; images.mkdir(exist_ok=True)
    current_figures = {path.name for path in (comparison/'figuras').glob('*')
                       if path.suffix in ('.png', '.svg')}
    for path in images.glob('resultado-*'):
        if path.name.removeprefix('resultado-') not in current_figures:
            archive = output/'historico'/'figuras-anteriores'
            archive.mkdir(parents=True, exist_ok=True)
            path.rename(archive/path.name)
    for path in (original/'imagens').glob('*'):
        if path.suffix in ('.png','.svg'): shutil.copy2(path,images/path.name)
    for path in (comparison/'figuras').glob('*'):
        if path.suffix in ('.png','.svg'): shutil.copy2(path,images/('resultado-'+path.name))
    for name in ('resultados-discussao.typ','conclusao.typ','analise-resultados.md','conclusao.md','pendencias-metodologicas.md'):
        shutil.copy2(comparison/name,output/name)
    report = (output/'analise-resultados.md').read_text()
    report = report.replace('(figuras/', '(imagens/resultado-').replace('(../manifest.json)', '(manifest.json)')
    report = report.replace('`../historico/avaliacao-m128`',
                            '`../../results/tcc/historico/avaliacao-m128`')
    archived = output/'resumos'
    archived.mkdir(exist_ok=True)
    for name in runs:
        shutil.copy2(root/name/'summary.json', archived/f'{name}-summary.json')
        report = report.replace(f'(../{name}/summary.json)', f'(resumos/{name}-summary.json)')
    (output/'analise-resultados.md').write_text(report)
    shutil.copy2(root/'manifest.json', output/'manifest.json')
    shutil.copytree(root/'source', output/'source', dirs_exist_ok=True)
    if (root/'statistical-analysis-source').exists():
        shutil.copytree(root/'statistical-analysis-source', output/'statistical-analysis-source', dirs_exist_ok=True)
    for name in ('funcionamento-codigo.md','testes-executados.txt','verificacao-integridade.json','fontes-e-implementacao.md','verificacao-entrega.md'):
        shutil.copy2(comparison/name, output/name)
    for name in ('auditoria.log', 'cli-help.txt', 'reavaliacao-nist.log'):
        if (comparison/name).exists():
            shutil.copy2(comparison/name, output/name)
    if (comparison/'reavaliacao-nist.json').exists():
        shutil.copy2(comparison/'reavaliacao-nist.json', output/'reavaliacao-nist.json')
    shutil.copy2(original/'referencias.bib',output/'referencias.bib')
    tables={'parameters':[],'energies':[],'times':[],'quantiles':[],'bits':[],'tests':[],
            'waiting_main':[], 'nist_main':[]}
    testnames={'monobit':'Frequência global','block_frequency':'Frequência em blocos','runs':'Runs'}
    states={'pass':'Sem reprovação','fail':'Reprovação','not_applicable':'Não aplicável'}
    for name,s in runs.items():
        params=s['parameters']
        tables['parameters'].append([name,pt(params['sigma_s'],3),pt(params['sigma_w'],3),pt(params['untruncated_sigma_E_eV'],6)])
        for title,key in (('Prévia','preflight'),('Oferecidas','realized_offered_barriers'),('Selecionadas','selected_barriers')):
            b=s[key]
            tables['energies'].append([name,title,f"{b['count']:,}".replace(',','.'),pt(b['mean_eV'],6),pt(b['std_eV'],6)])
        w=s['waiting_times']
        tables['waiting_main'].append([name,pt(w['mean_s']*1e12,3),pt(w['median_s']*1e12,3),pt(w['std_s']*1e12,3)])
        tables['times'].append([name,pt(w['mean_s']*1e12,6),pt(w['median_s']*1e12,6),pt(w['std_s']*1e12,6),pt(s['total_time_s']*1e6,6)])
        tables['quantiles'].append([name,*[pt(w['quantiles_s'][f'{q:02d}%']*1e12,5) for q in (1,5,25,75,95,99)]])
        tests=s['nist_sp800_22_subset']
        tables['nist_main'].append([name,*[p_value(tests['tests'][test]['p_value'])
                                         for test in ('monobit','block_frequency','runs')]])
        tables['bits'].append([name,pt(tests['ones_fraction'],3),str(s['digitization']['ties']),f"{tests['runs_count']:,}".replace(',','.'),f"{s['landscape']['unique_wells']:,}".replace(',','.'),f"{s['landscape']['unique_saddles']:,}".replace(',','.')])
        for test,t in tests['tests'].items():
            p=t['p_value']
            ptext='Não aplicável' if p is None else ('0 (numérico)' if p==0 else (pt(p,6) if p>.000001 else f'{p:.3e}'.replace('.',',')))
            tables['tests'].append([name,testnames[test],pt(t['statistic'],4) if t['statistic'] is not None else 'Não aplicável',ptext,states[t['status']]])
    data={'tables':tables,'raw_summaries':runs,'manifest':manifest}
    (output/'dados-resultados.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    source=(original/'tcc-final.typ').read_text()
    methodology = '= Metodologia\n\n'
    code_reference = ('As simulações deste trabalho foram executadas com o código Apeiron, '
                      'desenvolvido em Python e disponível no repositório '
                      '#link("https://github.com/IsaqueFranklin/apeiron-numbers")[`apeiron-numbers`] '
                      '@apeiron2026. A versão utilizada foi preservada com os dados, junto das sementes '
                      'e dos registros de execução.\n\n')
    if '@apeiron2026' not in source:
        if source.count(methodology) != 1:
            raise ValueError('methodology heading changed; inspect before adding the code reference')
        source = source.replace(methodology, methodology + code_reference, 1)
    # Present tense states the now implemented protocol; preserve personal fields and references.
    replacements={'Serão':'São','serão':'são','Será':'É','será':'é',
                  'será':'é','receberá':'recebe','produzirá':'produz','produzirão':'produzem',
                  'permitirão':'permitem','permanecerá':'permanece','continuarão':'continuam',
                  'terá':'tem','terão':'têm','definirá':'define','selecionará':'seleciona',
                  'atualizará':'atualiza','utilizará':'utiliza','investigará':'investiga',
                  'fornecerá':'fornece','se deslocará':'se desloca','distribuirá':'distribui',
                  'determinará':'determina','guardará':'guarda','trocarão':'trocam',
                  'dependerão':'dependem','apresentarão':'apresentam','ampliará':'amplia',
                  'poderá':'pode','poderão':'podem','mostrarão':'mostram',
                  'deverão':'devem','deverá':'deve','recuperará':'recupera',
                  'existirão':'existem','reconstruirá':'reconstrói'}
    for old,new in replacements.items(): source=source.replace(old,new)
    source = source.replace('Procedimentos previstos', 'Procedimentos adotados')
    obsolete='O @fluxo:metodologia reúne as etapas do procedimento. A implementação consultada já contém o modelo energético, o KMC, a conversão por mediana e os testes estatísticos. A execução conjunta de H0, H1 e H2 depende da extensão dos parâmetros de entrada, pois os desvios padrão ainda estão fixados no código. As análises são realizadas depois dessa adequação, sem antecipar resultados.'
    updated='O @fluxo:metodologia reúne as etapas do procedimento. A implementação foi estendida para receber uma configuração comum de parâmetros na paisagem e na verificação prévia, cadastrar H0, H1 e H2 em um comando e registrar a proveniência de cada rodada. As três execuções completas e as verificações científicas estão documentadas no capítulo de resultados.'
    if obsolete not in source: raise ValueError('source paragraph changed; inspect before replacing')
    source=source.replace(obsolete,updated)
    source=source.replace('Essas saídas comparativas são acrescentadas na adequação do programa.','Essas saídas comparativas foram acrescentadas ao programa e são apresentadas no capítulo de resultados.')
    paragraphs = source.split('\n\n')
    for i, paragraph in enumerate(paragraphs):
        if paragraph.startswith('Na comparação descritiva são examinados'):
            paragraphs[i] = 'A comparação principal utiliza os tempos médios, as medianas e os desvios padrão das esperas, junto dos resultados dos três testes estatísticos. Os gráficos apresentam um mesmo trecho das séries temporais, o equilíbrio por blocos, as contagens de runs e os valores-p. Quantis e estatísticas detalhadas das barreiras permanecem nos arquivos de apoio, sem ampliar a discussão principal.'
        elif paragraph.startswith('Para um milhão de saltos, cada rodada oferecerá'):
            paragraphs[i] = 'Os registros completos preservam doze milhões de barreiras oferecidas e um milhão de barreiras selecionadas por cenário. Esses conjuntos são mantidos separadamente da verificação prévia independente, pois a trajetória revisita ambientes e escolhe saltos conforme suas taxas.'
        elif paragraph.startswith('O código divide a sequência em blocos não sobrepostos'):
            paragraphs[i] = 'O código divide a sequência em blocos não sobrepostos de $M = 20000$ bits para $n = 10^6$. O tamanho automático é $M = max(20, ceil(n/50))$, buscando no máximo 50 blocos. Se $B = floor(n/M)$ e $pi_l$ é a proporção de uns no bloco $l$, a estatística implementada é'
        elif paragraph.startswith('Para $n = 10^6$, são usados'):
            paragraphs[i] = 'Para $n = 10^6$, são usados 50 blocos completos e nenhum bit residual. Todos os bits participam do teste. A rotina pode detectar trechos com concentração de símbolos mesmo quando o equilíbrio global é imposto pela mediana.'
        elif paragraph.startswith('O tamanho de $128$ bits corresponde'):
            paragraphs[i] = 'A configuração de 20000 bits atende às recomendações da seção 2.2.7 da SP 800-22. O tamanho é pelo menos 20 bits, supera 1% da sequência e resulta em menos de 100 blocos. A adequação desse parâmetro não representa conformidade com a bateria completa, pois somente três testes são utilizados.'
        elif paragraph.startswith('Os resultados são apresentados por cenário, preservando'):
            paragraphs[i] = 'Os resultados são apresentados por cenário, com parâmetros energéticos, tempos de espera, tamanho dos blocos, proporção de uns, número de runs, estatísticas, valores-p e condições de aplicação. A discussão principal relaciona os padrões binários ao comportamento temporal da vacância. Maior variação de energia não é tomada como garantia de uma sequência mais aleatória.'
    source = '\n\n'.join(paragraphs)
    source = source.replace('separator: [ – ]', 'separator: [ ]')
    # Author-provided full abstracts are written only after reading all complete summaries.
    abstract=(comparison/'resumos-atualizados.json')
    if abstract.exists():
        abstracts=json.loads(abstract.read_text())
        for start,key in (('#resumo[','pt'),('#resumo-en[','en')):
            begin=source.index(start)+len(start)
            end=source.index('\n]',begin)
            source=source[:begin]+'\n  '+abstracts[key]+'\n'+source[end:]
    bibliography='#abnt-bibliography(read("referencias.bib"), titulo: "REFERÊNCIAS")'
    source=source.replace(bibliography,'#include "resultados-discussao.typ"\n\n#include "conclusao.typ"\n\n'+bibliography)
    (output/'tcc-com-resultados.typ').write_text(source)
    standalone='''#import "@preview/abntyp:0.1.5": *
#show: normas-abnt.with(fonte: sys.inputs.at("fonte", default: "Times New Roman"))
#set math.equation(numbering: "(1)", number-align: right + horizon)
#set figure.caption(position: top, separator: [ ])
#show figure.caption: set text(size: 10pt)
#show figure: set align(center)
#textual()
#include "resultados-discussao.typ"
#include "conclusao.typ"
#abnt-bibliography(read("referencias.bib"), titulo: "REFERÊNCIAS")
'''
    (output/'resultados-conclusao.typ').write_text(standalone)
    (output/'README.md').write_text('''# Resultados, discussão e conclusão do TCC

Os textos foram escritos depois da leitura das três rodadas completas, sem geração automática de conclusões. As tabelas em `dados-resultados.json` são formatadas diretamente a partir de seus resumos; os gráficos vêm das amostras completas.

- `tcc-com-resultados.typ`: versão unificada, incluindo o material anterior, a metodologia atualizada para a implementação e os capítulos novos. Os originais foram preservados.
- `resultados-conclusao.typ`: documento de conferência com os capítulos novos e referências.
- `resultados-discussao.typ` e `conclusao.typ`: fragmentos para integração.
- `analise-resultados.md` e `conclusao.md`: conteúdo em Markdown.
- `dados-resultados.json`: tabelas e resumos originais utilizados.
- `imagens`: gráficos dos resultados em PNG 350 dpi e SVG, mais as ilustrações anteriores necessárias ao TCC unificado.
- `referencias.bib`: bibliografia existente, conservada.
- `funcionamento-codigo.md`: funcionamento, CLI, sementes, saídas e limites de interpretação.
- `resumos`, `manifest.json` e `source`: resumos completos, proveniência e código efetivamente executado.
- `testes-executados.txt`, `verificacao-integridade.json` e `verificacao-entrega.md`: testes e auditoria.

Na plataforma Typst, envie todos esses arquivos, preservando a pasta `imagens`, e compile `tcc-com-resultados.typ`. A biblioteca permanece `@preview/abntyp:0.1.5`. O padrão solicita Times New Roman; a conferência local usa Liberation Serif pelo parâmetro `fonte`, por disponibilidade.

A avaliação usa 50 blocos de 20.000 bits por cenário e os três testes da SP 800-22. A revisão estatística foi feita sobre os mesmos bits, sem repetir a simulação. A mediana da própria rodada e uma execução por cenário são limitações descritas em `pendencias-metodologicas.md`. Campos pessoais de aprovação e agradecimentos não foram preenchidos.

O capítulo principal contém três tabelas e quatro gráficos. A comparação de quantis e de barreiras oferecidas e selecionadas foi retirada do texto. Os registros completos e a avaliação anterior com 128 bits por bloco permanecem em `../../results/tcc`, incluindo seu histórico.

Os registros completos de eventos não integram este pacote de texto. Consulte `../../results/tcc/LEIA-ME.md` e o campo `delivery.all_event_records_persistent` do manifesto para conferir o armazenamento dos CSVs.
''')
    return output


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('results',type=Path)
    parser.add_argument('--thesis-dir', type=Path, default=THESIS_ROOT,
                        help='TCC source directory (default: ../tcc-apeiron-words)')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    output = args.output if args.output is not None else args.thesis_dir/'gerados'/'experimentos-tcc'
    print(assemble(args.results, output, args.thesis_dir))
