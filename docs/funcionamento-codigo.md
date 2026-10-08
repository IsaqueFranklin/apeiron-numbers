# Funcionamento do Apeiron e reprodução do protocolo do TCC

## Finalidade e módulos

O Apeiron simula uma única vacância em uma rede FCC periódica. A aproximação estatística associa energias a assinaturas químicas locais; não calcula energias com NEB ou um potencial interatômico e não prepara ordenamento de curto alcance. A saída científica principal é a sequência de tempos de espera. Os bits são uma transformação exploratória dessa sequência, reproduzível pelas sementes.

- `config.py`: configuração imutável e compartilhada de médias e desvios dos poços/selas; validação e rejeição finita das selas.
- `model.py`: rede, ocupação química, chaves de energia, caches, taxas, escolha do salto e relógio KMC.
- `analysis.py`: verificação prévia independente, momentos em lotes, extração pela mediana e os três testes.
- `cli.py`: entrada, cenários independentes, escrita incremental, proveniência e estados da execução.
- `plotting.py`: diagnóstico individual, conservando a opção de replotar dados existentes.
- `comparison.py`: tabelas e figuras comparativas, incluindo somente rodadas concluídas.

## Rede, espécies e configuração inicial

Cada célula convencional possui quatro sítios FCC. Com oito células por eixo, existem 2.048 sítios, 2.047 átomos e uma vacância. A rede é representada por coordenadas inteiras com soma par num cubo de lado 16. As doze direções de salto são permutações de `(±1, ±1, 0)`.

O mapa `occupancy` associa cada coordenada a um código: Co=0, Ni=1, Cr=2, Fe=3, Mn=4, vacância=5. Um embaralhamento com a semente da composição distribui 410 Co, 410 Ni e 409 de cada outro elemento. A vacância começa em `(0,0,0)`. A cada salto, somente as duas ocupações envolvidas são alteradas. As espécies são conservadas.

## Energias e caches

Uma configuração comum `EnergyConfig` contém `mu_s`, `mu_w`, `sigma_s` e `sigma_w`. O mesmo objeto é passado à paisagem e à verificação prévia. Os poços têm distribuição nominal N(mu_w,sigma_w²); as selas têm N(mu_s,sigma_s²) antes da rejeição.

A chave do poço é a sequência ordenada dos códigos dos doze vizinhos, sem coordenadas absolutas. A chave da sela contém a espécie móvel, a direção canônica da ligação e os códigos dos 18 sítios externos da união das duas primeiras vizinhanças. Os dois extremos são omitidos dessa vizinhança; o átomo móvel é registrado separadamente. A orientação canônica é escolhida comparando as duas descrições geométricas da mesma ponte, o que garante a mesma chave na reversão imediata.

As assinaturas exatas usam bytes para reduzir memória. Antes de obter a semente local por BLAKE2b, o código restaura a representação original da tupla do commit de referência. Essa compactação conserva os sorteios daquele código: a comparação independente de 100 eventos de H2 é um teste automatizado.

Ao encontrar uma chave nova, o programa obtém uma energia reproduzível pela chave e pela semente e a coloca no cache. Visitas posteriores reutilizam a energia. Uma mesma configuração em outra coordenada tem a mesma energia; uma configuração rotacionada geralmente tem outra chave. O código preserva a simetria translacional e não agrupa rotações.

Para uma sela nova, são calculados o poço atual e o poço que existiria no destino depois da troca virtual. A sela é aceita somente se `s >= max(w_origem,w_destino)`. Apenas candidatos energéticos são rejeitados; todo evento KMC selecionado é executado. O código registra candidatos e rejeições e limita a 100.000 candidatos por sorteio para interromper parâmetros extremos em vez de entrar em um laço infinito.

Sigma zero produz exatamente a média. Uma sela constante com poços Gaussianos de desvio positivo é recusada, pois os poços possuem suporte não limitado superiormente; a aceitação não seria possível para todas as configurações. Uma sela constante abaixo de um poço constante também é recusada. Médias e sigmas devem ser finitos, e sigmas não podem ser negativos.

## Passo KMC e relógio

1. Obter o poço da configuração atual, compartilhado pelas doze saídas.
2. Para cada vizinho, obter o poço de destino após a troca virtual e a sela reversível.
3. Calcular `E_ij = s_ij - w_i` e `r_ij = nu*exp(-E_ij/(k_B*T))`, com `k_B=8.617333262e-5 eV/K`.
4. Somar as taxas com `math.fsum` para obter R da configuração anterior ao salto.
5. Sortear `u1` para escolher o destino pela soma acumulada das taxas.
6. Sortear outro uniforme independente `u2` em (0,1] e calcular `dt=-log(u2)/R`.
7. Acrescentar dt ao relógio e efetuar a troca da vacância com o átomo.

O tempo de espera é sempre um novo sorteio, mesmo quando as energias são reutilizadas. Em H0 as taxas permanecem constantes e os destinos equiprováveis, mas os tempos continuam estocásticos. Para os parâmetros solicitados, `R0≈7.45435e10 s⁻¹` e `1/R0≈1.34150e-11 s`.

## Cenários e comando

```bash
python main.py \
  --scenarios H0:0:0 H1:0.152:0.055 H2:0.304:0.110 \
  --jumps 1000000 --temperature 1273 --cells 8 \
  --attempt-frequency 1e13 --seed 2020 \
  --preflight-samples 1000000 --output results/tcc --plot
```

A ordem no cadastro é **nome:sigma_s:sigma_w**. H1 é metade dos desvios de H2, uma escolha do planejamento. H2 conserva os valores arredondados derivados de Thomas e Patala. A largura nominal da barreira antes da restrição é `hypot(sigma_s,sigma_w)`, aproximadamente 0,323289344 eV em H2; ela não é ajustada para 0,320.

No modo de uma rodada, use `--sigma-s 0.304 --sigma-w 0.110`, sem `--scenarios`. As médias podem ser configuradas por `--mu-s` e `--mu-w`; os padrões são 0,81 e 0 eV. Os demais padrões estão em `python main.py --help`. Nomes de cenário devem conter letras, números, hífen ou sublinhado; nomes repetidos e caminhos são recusados.

Cada cenário inicializa uma nova rede e novos caches. As sementes são base para composição, base+1 para energias, base+2 para KMC e base+10 para a verificação independente. Compartilhar a semente química controla a condição inicial; as probabilidades diferentes produzem trajetórias diferentes, sem forçar os mesmos saltos.

## As três amostras energéticas

A verificação prévia sorteia um milhão de pares de poços independentes e uma sela aceita por par, com a mesma configuração do cenário. Os valores estão em `preflight_barriers.npy` e as estatísticas em `summary.json/preflight`.

As doze barreiras oferecidas por evento totalizam doze milhões por rodada; são correlacionadas pela paisagem visitada. Os valores são registrados em `events.csv`, com momentos em `realized_offered_barriers`. A barreira selecionada, identificada pelo índice `selected_neighbor` de base zero, totaliza um milhão; sua estatística está em `selected_barriers`. Sua distribuição é ponderada pela escolha das taxas e não substitui o desvio nominal do modelo.

O KMC agrega momentos populacionais e histogramas em lotes de 4.096 eventos. Mantém um array NumPy com um milhão de tempos para a mediana e os quantis, sem manter doze milhões de objetos float em listas. Os histogramas completos por bins comuns ficam em `barrier_histograms.npz`; os valores individuais permanecem no registro incremental.

## Bits e testes

Para cada rodada, `theta=median(dt)` e `bit=1` somente se `dt>theta`; empates produzem zero. A ordem dos eventos é preservada. O limiar e o número de empates são registrados. A quebra de linha final de `bits.txt` não é um bit.

O código executa frequência global, frequência em blocos e runs da **SP 800-22 Rev. 1a**, com alpha=0,01. Cada resultado guarda estatística, valor-p e estado `pass`, `fail` ou `not_applicable`; valores ausentes são `null`. Valores numericamente iguais a zero continuam sendo zero. Entradas com menos de 100 bits são insuficientes. O runs exige `abs(proporcao_uns-0.5)<2/sqrt(n)` e guarda os dois lados dessa condição.

O tamanho automático dos blocos é `M=max(20,ceil(n/50))`. Com um milhão de bits, são 50 blocos de 20.000 bits e nenhum resíduo. A configuração atende a M≥20, M>0,01*n e B<100 da seção 2.2.7 da SP 800-22. `--block-size` permite definir M explicitamente, com validação dessas recomendações pela CLI para sequências com pelo menos 100 bits. Os resíduos, quando houver, ficam fora somente do teste em blocos. A mediana dos próprios dados força o equilíbrio global quando não há empates, de modo que o monobit não fornece evidência independente de qualidade. A aprovação em três testes não é aprovação na bateria NIST completa, avaliação SP 800-90B ou certificação criptográfica.

## Arquivos, rastreabilidade e falhas

`manifest.json` registra comando, sementes, cenários, duração, versões, commit e hashes dos arquivos executados. A pasta `source` conserva uma cópia do código efetivamente usado; isso é necessário porque o commit pode ter alterações locais. `summary.json` guarda composição inicial/final, número de chaves, rejeições, tempos, quantis e resultados dos testes.

A gravação dos CSVs é incremental. Resumos e manifesto são atualizados atomicamente. Antes de executar, o código verifica espaço; se uma rodada falha, guarda `incomplete`, o erro e o número de eventos executados. A comparação exclui rodadas incompletas. Arquivos de rodada existentes não são sobrescritos; uma nova execução deve usar outro destino.

`--plot` produz PNG de 350 dpi, SVG e o diagnóstico PDF individual. Os quatro gráficos principais mostram o mesmo trecho dos primeiros 2.000 saltos, as frações de uns nos blocos usados no teste, as contagens de runs e os valores-p. As escalas são comuns entre cenários. Essa redução da série temporal não reduz a amostra dos testes. Zeros são contados e omitidos apenas dos eixos logarítmicos. Os valores-p usam eixo linear e conservam zeros numéricos. Quantis e estatísticas das barreiras ficam nos arquivos de apoio.

Para replotar uma rodada: `python -m apeiron.plotting DIRETORIO`. Para recompor a comparação: `python -c 'from apeiron.comparison import compare_results; compare_results("results/tcc")'`.

## Verificação e leitura científica

```bash
python -m unittest discover -s tests -v
```

Os testes cobrem geometria, composição, reversão em todas as direções, igualdade por translação, propagação de parâmetros, limite homogêneo, preservação da trajetória de referência, caches independentes, reprodução das três rodadas, empates, aplicabilidade e correspondência evento–tempo–bit. As verificações pequenas não substituem as três execuções completas do protocolo.

Uma realização por cenário permite uma comparação descritiva. A alteração simultânea dos dois desvios não permite atribuir efeitos isolados a poços ou selas. O tempo necessário para um número fixo de saltos não é um coeficiente de difusão: a análise de transporte exigiria deslocamentos e uma escala física de rede. As conclusões tratam do modelo numérico, sem representar medições de uma liga física.

Fontes: [Thomas e Patala, manuscrito público](https://arxiv.org/pdf/2003.06900), [artigo](https://doi.org/10.1016/j.actamat.2020.06.022) e [SP 800-22 Rev. 1a](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf). O índice e o truncamento do Apeiron são adaptações do artigo: o código conserva suas assinaturas exatas e a restrição local de sela acima dos dois poços, em vez de reproduzir literalmente todos os detalhes de indexação e truncamento do suplemento S1.

## Reavaliação estatística de resultados existentes

As três simulações completas foram mantidas. `python scripts/reanalyse_nist.py results/tcc` recalculou os testes sobre os mesmos bits e gerou os gráficos atuais. A avaliação anterior foi preservada em `results/tcc/historico/avaliacao-m128`. O script recusa sobrescrever esse histórico.

`comparacao/reavaliacao-nist.json` registra o protocolo novo, comando, duração, versões, hashes do código e SHA-256 dos bits e tempos antes e depois. `statistical-analysis-source` contém o código da reavaliação. `source` continua sendo o código efetivamente usado na simulação original. Os hashes confirmam que os bits, tempos e eventos não foram alterados. A auditoria pode ser reproduzida com `python scripts/verify_tcc_results.py results/tcc`.

## Organização local dos arquivos

O repositório `random-heas` contém o código, os testes, os scripts e a documentação. O material do TCC foi separado na pasta irmã `../tcc-apeiron-words`, com os documentos finais em `gerados/experimentos-tcc` e todos os dados em `results`.

Os comandos anteriores com `results/tcc` registram a localização usada nas execuções originais. A localização atual dos dados é `../tcc-apeiron-words/results/tcc`. Para conferir os mesmos resultados a partir da raiz do código:

```bash
python scripts/verify_tcc_results.py ../tcc-apeiron-words/results/tcc
```

Para novas rodadas, `--output` aceita um caminho externo. O script `build_thesis_bundle.py` usa `../tcc-apeiron-words` como diretório padrão do material do TCC e aceita `--thesis-dir` para outra localização. O simulador e seus testes funcionam sem essa pasta.
