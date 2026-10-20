# Résultats

Calibrage sur ARL0 = 370,4 (Shewhart 3σ) : EWMA λ = 0,2 → L = 2,859 ; CUSUM k = 0,5 → h = 4,755. 4 000 séries simulées par point, sous-groupes de taille 4.

## ARL selon le décalage de la moyenne (intervalle de confiance à 95 %)

|      | Shewhart exact   | Shewhart 3σ   | Shewhart + règles WE   | EWMA (λ = 0,2)   | CUSUM (k = 0,5)   |
|:-----|:-----------------|:--------------|:-----------------------|:-----------------|:------------------|
| 0,00 | 370,4            | 375,4 ± 11,6  | 91,4 ± 2,6             | 360,7 ± 11,2     | 355,3 ± 10,7      |
| 0,25 | 155,2            | 156,1 ± 4,7   | 26,6 ± 0,7             | 36,3 ± 1,0       | 34,9 ± 0,9        |
| 0,50 | 43,9             | 43,8 ± 1,4    | 9,1 ± 0,2              | 9,7 ± 0,2        | 9,8 ± 0,2         |
| 0,75 | 15,0             | 15,0 ± 0,4    | 4,8 ± 0,1              | 5,2 ± 0,1        | 5,5 ± 0,1         |
| 1,00 | 6,3              | 6,3 ± 0,2     | 3,1 ± 0,1              | 3,6 ± 0,0        | 3,9 ± 0,0         |
| 1,50 | 2,0              | 2,0 ± 0,0     | 1,7 ± 0,0              | 2,3 ± 0,0        | 2,5 ± 0,0         |
| 2,00 | 1,2              | 1,2 ± 0,0     | 1,2 ± 0,0              | 1,8 ± 0,0        | 1,9 ± 0,0         |

## Décalage de +0,5σ : délais de détection sur 4 000 simulations

| méthode              |   délai médian | détecté en 5 sous-groupes   | en 10   | en 20   |
|:---------------------|---------------:|:----------------------------|:--------|:--------|
| Shewhart 3σ          |             31 | 11 %                        | 20 %    | 37 %    |
| Shewhart + règles WE |              8 | 32 %                        | 71 %    | 94 %    |
| EWMA (λ = 0,2)       |              8 | 23 %                        | 66 %    | 94 %    |
| CUSUM (k = 0,5)      |              9 | 18 %                        | 65 %    | 95 %    |

## Scénario simulé : délai de détection (sous-groupes après le début du défaut)

| caractéristique             | défaut                           |   Shewhart 3σ |   Shewhart + règles WE |   EWMA (λ = 0,2) |   CUSUM (k = 0,5) | carte R   |
|:----------------------------|:---------------------------------|--------------:|-----------------------:|-----------------:|------------------:|:----------|
| Humidité                    | décalage de +0,5σ                |             2 |                      2 |                4 |                 4 | —         |
| Résistance à la compression | dérive de −0,05σ par sous-groupe |            20 |                     15 |               15 |                20 | —         |
| Granulométrie               | écart-type × 2                   |             1 |                      1 |               28 |                 2 | 5         |

Fausses alarmes avant le défaut : Humidité : Shewhart 3σ 1, Shewhart + règles WE 1, EWMA (λ = 0,2) 0, CUSUM (k = 0,5) 0; Résistance à la compression : Shewhart 3σ 1, Shewhart + règles WE 2, EWMA (λ = 0,2) 1, CUSUM (k = 0,5) 0; Granulométrie : Shewhart 3σ 0, Shewhart + règles WE 1, EWMA (λ = 0,2) 0, CUSUM (k = 0,5) 0.

## Capabilité (phase I)

| caractéristique             | tolérance              | σ estimé (R̄/d2)   | Cp   | Cpk   |
|:----------------------------|:-----------------------|:-------------------|:-----|:------|
| Humidité                    | 3,6 – 5,4 %            | 0,24               | 1,26 | 1,20  |
| Résistance à la compression | 90,0 – 150,0 kPa       | 8,28               | 1,21 | 1,17  |
| Granulométrie               | 43,0 – 67,0 indice AFS | 4,07               | 0,98 | 0,97  |
