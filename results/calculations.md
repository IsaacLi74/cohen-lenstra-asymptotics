# Residual-ratio estimates

The final dyadic checkpoint is $X=2^{30}$ (1,073,741,824), with 326,378,404 fields.
The [README](../README.md#residual-ratio-method) gives the fixed main terms and definitions.
The [input manifest](../data/processed/p3_statistics.manifest.json) records the raw shards and SHA-256 hashes.

Each exponent column uses the endpoint and one earlier dyadic checkpoint. A dash marks a missing checkpoint, a zero residual, or opposite signs. The gap 4 coefficient is $\widehat b_{n-4,n}=E(2^n)/2^{n\widehat\theta_{n-4,n}}$. Each coefficient is paired with its own effective exponent. The $b_{5/6}=E(X)/X^{5/6}$ comparison is reported for $p=3$.
Types with fewer than 10,000 observations at the endpoint are marked **Sparse**. Within each prime, the Sylow types are ordered by $\lvert\mathrm{Aut}(G)\rvert$.

| Statistic | Count or sum | Observed / $N(X)$ | Leading mean or proportion | $\widehat\theta_{29,30}$ | $\widehat\theta_{28,30}$ | $\widehat\theta_{26,30}$ | $\widehat b_{26,30}$ |
|---|---:|---:|---:|---:|---:|---:|---:|
| [$3\mid h_D$](../figures/statistics/divisible.png) | 139,775,918 | 0.428263379 | 0.439873922 | 0.82939 | 0.82850022 | 0.82756928 | -0.12731433 |
| [3-torsion](../figures/statistics/torsion.png) | 638,841,092 | 1.957363245 | 2.000000000 | 0.8323687 | 0.83333374 | 0.83330132 | -0.41499525 |
| [$5\mid h_D$](../figures/statistics/p5/divisible.png) | 77,719,065 | 0.238125636 | 0.239667204 | 0.73242938 | 0.75463528 | 0.73253653 | -0.12195922 |
| [5-torsion](../figures/statistics/p5/torsion.png) | 650,168,024 | 1.992068152 | 2.000000000 | 0.72945526 | 0.74484817 | 0.72838474 | -0.68409638 |
| [$7\mid h_D$](../figures/statistics/p7/divisible.png) | 53,062,311 | 0.162579112 | 0.163204593 | 0.76133332 | 0.69941972 | 0.70283706 | -0.091764651 |
| [7-torsion](../figures/statistics/p7/torsion.png) | 651,267,184 | 1.995435899 | 2.000000000 | 0.74479821 | 0.69637599 | 0.68886249 | -0.89539185 |

## Sylow 3

46 types in the data. [Compare the common groups](../figures/statistics/sylow_comparison.png).

| Group | Aut size | Count | Observed / $N(X)$ | CL probability | Gap 1 exponent | Gap 2 exponent | Gap 4 exponent | Gap 4 coefficient | $b_{5/6}$ | Status |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| [trivial](../figures/statistics/sylow_trivial.png) | 1 | 186,602,486 | 0.57173662 | 0.56012608 | 0.82938165 | 0.82848148 | 0.82760258 | 0.12722711 | 0.11293427 |  |
| [C3](../figures/statistics/sylow_C3.png) | 2 | 89,635,292 | 0.2746361 | 0.28006304 | 0.83658955 | 0.8422592 | 0.83799694 | -0.047908002 | -0.052786716 |  |
| [C9](../figures/statistics/sylow_C9.png) | 6 | 29,858,194 | 0.091483363 | 0.093354346 | 0.8354881 | 0.82761965 | 0.83359692 | -0.018099194 | -0.018198672 |  |
| [C27](../figures/statistics/sylow_C27.png) | 18 | 9,952,379 | 0.030493375 | 0.031118115 | 0.8620148 | 0.81005686 | 0.8353365 | -0.0058287999 | -0.0060767243 |  |
| [C3xC3](../figures/statistics/sylow_C3xC3.png) | 48 | 3,232,964 | 0.0099055696 | 0.011669293 | 0.83059444 | 0.8399611 | 0.84141855 | -0.014500571 | -0.017155439 |  |
| [C81](../figures/statistics/sylow_C81.png) | 54 | 3,317,214 | 0.010163706 | 0.010372705 | 0.8502674 | 0.8585453 | 0.83266178 | -0.002061484 | -0.0020328962 |  |
| [C9xC3](../figures/statistics/sylow_C9xC3.png) | 108 | 1,430,639 | 0.0043833752 | 0.0051863526 | 0.84142436 | 0.84619601 | 0.8489891 | -0.0056401277 | -0.0078104239 |  |
| [C243](../figures/statistics/sylow_C243.png) | 162 | 1,104,084 | 0.0033828341 | 0.0034575684 | 0.71526056 | 0.63895545 | 0.54849475 | -0.27154286 | -0.00072692511 |  |
| [C27xC3](../figures/statistics/sylow_C27xC3.png) | 324 | 478,861 | 0.0014671957 | 0.0017287842 | 0.85720144 | 0.83324095 | 0.83828899 | -0.0022952807 | -0.0025444263 |  |
| [C729](../figures/statistics/sylow_C729.png) | 486 | 354,950 | 0.0010875413 | 0.0011525228 | 0.20413394 | 0.28807569 | 0.36651169 | -10.389149 | -0.00063206314 |  |
| [C81xC3](../figures/statistics/sylow_C81xC3.png) | 972 | 159,457 | 0.0004885648 | 0.0005762614 | 0.82340488 | 0.81054411 | 0.79479056 | -0.0019012094 | -0.00085300982 |  |
| [C2187](../figures/statistics/sylow_C2187.png) | 1,458 | 90,877 | 0.0002784406 | 0.00038417426 | 0.54196584 | 0.57774952 | 0.65383849 | -0.042971906 | -0.0010284535 |  |
| [C243xC3](../figures/statistics/sylow_C243xC3.png) | 2,916 | 51,675 | 0.00015832849 | 0.00019208713 | 0.64898667 | 0.65168849 | 0.64292286 | -0.017216079 | -0.00032836457 |  |
| [C9xC9](../figures/statistics/sylow_C9xC9.png) | 3,888 | 39,385 | 0.00012067281 | 0.00014406535 | 0.90476415 | 0.86926765 | 0.82057736 | -0.00029665152 | -0.00022753522 |  |
| [C6561](../figures/statistics/sylow_C6561.png) | 4,374 | 14,719 | 4.5097959e-05 | 0.00012805809 | 0.80312113 | 0.82281091 | 0.86636708 | -0.00040599101 | -0.00080693934 |  |
| [C27xC9](../figures/statistics/sylow_C27xC9.png) | 8,748 | 17,554 | 5.3784196e-05 | 6.4029044e-05 | 0.91050044 | 0.79516649 | 0.8017399 | -0.00019221795 | -9.9649892e-05 |  |
| [C729xC3](../figures/statistics/sylow_C729xC3.png) | 8,748 | 13,229 | 4.0532706e-05 | 6.4029044e-05 | 0.65184944 | 0.6742677 | 0.73111863 | -0.0019145306 | -0.00022854494 |  |
| [C3xC3xC3](../figures/statistics/sylow_C3xC3xC3.png) | 11,232 | 7,345 | 2.2504553e-05 | 4.9868775e-05 | 0.88143317 | 0.88545496 | 0.89527361 | -7.3413869e-05 | -0.00026616722 | Sparse |
| [C19683](../figures/statistics/sylow_C19683.png) | 13,122 | 1,320 | 4.0443852e-06 | 4.2686029e-05 | 0.93403953 | 0.94538224 | 0.96409812 | -2.4780117e-05 | -0.00037586084 | Sparse |
| [C9xC3xC3](../figures/statistics/sylow_C9xC3xC3.png) | 23,328 | 3,492 | 1.0699237e-05 | 2.4010892e-05 | 0.88883675 | 0.88597815 | 0.87147218 | -5.858351e-05 | -0.00012948024 | Sparse |
| [C2187xC3](../figures/statistics/sylow_C2187xC3.png) | 26,244 | 2,202 | 6.7467699e-06 | 2.1343015e-05 | 0.85327287 | 0.8440166 | 0.87905233 | -5.4869243e-05 | -0.00014197524 | Sparse |
| [C81xC9](../figures/statistics/sylow_C81xC9.png) | 26,244 | 5,651 | 1.7314258e-05 | 2.1343015e-05 | 0.69371424 | 0.66991313 | 0.61592125 | -0.0036022106 | -3.9187029e-05 | Sparse |
| [C59049](../figures/statistics/sylow_C59049.png) | 39,366 | 6 | 1.8383569e-08 | 1.4228676e-05 | 0.99813482 | 0.99906741 | 0.99953371 | -4.3614972e-06 | -0.00013822115 | Sparse |
| [C27xC3xC3](../figures/statistics/sylow_C27xC3xC3.png) | 69,984 | 1,155 | 3.5388371e-06 | 8.0036305e-06 | 0.89409017 | 0.88367472 | 0.90578426 | -9.6266633e-06 | -4.34283e-05 | Sparse |
| [C243xC9](../figures/statistics/sylow_C243xC9.png) | 78,732 | 1,438 | 4.4059288e-06 | 7.1143382e-06 | 0.64578348 | 0.66525059 | 0.7479588 | -0.0001554872 | -2.6344245e-05 | Sparse |
| [C6561xC3](../figures/statistics/sylow_C6561xC3.png) | 78,732 | 191 | 5.8521029e-07 | 7.1143382e-06 | 0.95402635 | 0.94432053 | 0.96904018 | -3.7780846e-06 | -6.3507741e-05 | Sparse |
| [C81xC3xC3](../figures/statistics/sylow_C81xC3xC3.png) | 209,952 | 408 | 1.2500827e-06 | 2.6678768e-06 | 0.82963765 | 0.87441454 | 0.83718986 | -1.2727898e-05 | -1.3790647e-05 | Sparse |
| [C729xC9](../figures/statistics/sylow_C729xC9.png) | 236,196 | 240 | 7.3534277e-07 | 2.3714461e-06 | 0.82228766 | 0.84065576 | 0.88921553 | -4.9786963e-06 | -1.5914104e-05 | Sparse |
| [C27xC27](../figures/statistics/sylow_C27xC27.png) | 314,928 | 466 | 1.4277905e-06 | 1.7785846e-06 | 0.70475775 | 0.54046518 | 0.64250252 | -0.00018046686 | -3.4121139e-06 | Sparse |
| [C243xC3xC3](../figures/statistics/sylow_C243xC3xC3.png) | 629,856 | 101 | 3.0945675e-07 | 8.8929228e-07 | 0.91848929 | 0.89715368 | 0.86619792 | -2.8476033e-06 | -5.6399635e-06 | Sparse |
| [C9xC9xC3](../figures/statistics/sylow_C9xC9xC3.png) | 629,856 | 127 | 3.8911888e-07 | 8.8929228e-07 | 0.84116028 | 0.80388945 | 0.79244264 | -1.1386016e-05 | -4.8651031e-06 | Sparse |
| [C2187xC9](../figures/statistics/sylow_C2187xC9.png) | 708,588 | 14 | 4.2894995e-08 | 7.9048203e-07 | 0.97654067 | 0.98247626 | 0.97987722 | -3.4531074e-07 | -7.2716546e-06 | Sparse |
| [C81xC27](../figures/statistics/sylow_C81xC27.png) | 708,588 | 141 | 4.3201388e-07 | 7.9048203e-07 | 1.089008 | 0.86030472 | 0.92013274 | -5.7351566e-07 | -3.4867597e-06 | Sparse |
| [C27xC9xC3](../figures/statistics/sylow_C27xC9xC3.png) | 1,417,176 | 68 | 2.0834712e-07 | 3.9524101e-07 | 0.74090547 | 0.83197078 | 0.77763448 | -5.7886505e-06 | -1.8178857e-06 | Sparse |
| [C729xC3xC3](../figures/statistics/sylow_C729xC3xC3.png) | 1,889,568 | 25 | 7.6598205e-08 | 2.9643076e-07 | 0.90281023 | 0.81481257 | 0.89217746 | -6.2899624e-07 | -2.1382746e-06 | Sparse |
| [C243xC27](../figures/statistics/sylow_C243xC27.png) | 2,125,764 | 22 | 6.7406421e-08 | 2.6349401e-07 | 0.87071456 | 0.97780251 | 0.89343201 | -5.4660864e-07 | -1.9073113e-06 | Sparse |
| [C81xC9xC3](../figures/statistics/sylow_C81xC9xC3.png) | 4,251,528 | 19 | 5.8214636e-08 | 1.31747e-07 | 0.9411048 | 0.72783382 | 0.9575177 | -5.406963e-08 | -7.1523707e-07 | Sparse |
| [C2187xC3xC3](../figures/statistics/sylow_C2187xC3xC3.png) | 5,668,704 | 2 | 6.1278564e-09 | 9.8810253e-08 | 1 | 0.9538173 | 0.97690865 | -4.5535955e-08 | -9.0150625e-07 | Sparse |
| [C729xC27](../figures/statistics/sylow_C729xC27.png) | 6,377,292 | 2 | 6.1278564e-09 | 8.7831336e-08 | 0.89566176 | 0.94783088 | 0.97391544 | -4.2719789e-08 | -7.9471615e-07 | Sparse |
| [C243xC9xC3](../figures/statistics/sylow_C243xC9xC3.png) | 12,754,584 | 2 | 6.1278564e-09 | 4.3915668e-08 | 1.2552618 | 0.8915926 | 0.9457963 | -3.5455322e-08 | -3.6755575e-07 | Sparse |
| [C3xC3xC3xC3](../figures/statistics/sylow_C3xC3xC3xC3.png) | 24,261,120 | 1 | 3.0639282e-09 | 2.3087396e-08 | 0.79458637 | 0.89729318 | 0.94864659 | -1.7706319e-08 | -1.9476493e-07 | Sparse |
| [C81xC81](../figures/statistics/sylow_C81xC81.png) | 25,509,168 | 2 | 6.1278564e-09 | 2.1957834e-08 | 1 | 1.3531436 | 0.88198085 | -5.5991342e-08 | -1.5397555e-07 | Sparse |
| [C729xC9xC3](../figures/statistics/sylow_C729xC9xC3.png) | 38,263,752 | 1 | 3.0639282e-09 | 1.4638556e-08 | 1.4436162 | 0.83059626 | 0.91529813 | -2.0476826e-08 | -1.1258448e-07 | Sparse |
| [C27xC27xC3](../figures/statistics/sylow_C27xC27xC3.png) | 51,018,336 | 2 | 6.1278564e-09 | 1.0978917e-08 | -0.17836379 | 0.4108181 | 0.70540905 | -6.7465154e-07 | -4.7185455e-08 | Sparse |
| [C81xC27xC3](../figures/statistics/sylow_C81xC27xC3.png) | 114,791,256 | 1 | 3.0639282e-09 | 4.8795187e-09 | — | 0.28684989 | 0.64342494 | -9.1629005e-07 | -1.7659944e-08 | Sparse, Opposite signs for gap 1 |
| [C9xC9xC9](../figures/statistics/sylow_C9xC9xC9.png) | 221,079,456 | 1 | 3.0639282e-09 | 2.5335962e-09 | -1.7607267 | — | — | — | 5.1584533e-09 | Sparse, Opposite signs for gap 2, Opposite signs for gap 4 |

## Sylow 5

19 types in the data. [Compare the common groups](../figures/statistics/p5/sylow_comparison.png).

| Group | Aut size | Count | Observed / $N(X)$ | CL probability | Gap 1 exponent | Gap 2 exponent | Gap 4 exponent | Gap 4 coefficient | Status |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| [trivial](../figures/statistics/p5/sylow_trivial.png) | 1 | 248,659,339 | 0.76187436 | 0.7603328 | 0.7323755 | 0.75451153 | 0.73273362 | 0.12146678 |  |
| [C5](../figures/statistics/p5/sylow_C5.png) | 4 | 61,737,128 | 0.18915813 | 0.1900832 | 0.73819106 | 0.78471648 | 0.76784075 | -0.035123375 |  |
| [C25](../figures/statistics/p5/sylow_C25.png) | 20 | 12,352,741 | 0.037847912 | 0.03801664 | 0.96003003 | 0.93112704 | 0.9655159 | -0.00010505667 |  |
| [C125](../figures/statistics/p5/sylow_C125.png) | 100 | 2,468,429 | 0.0075630893 | 0.007603328 | 0.75257288 | 0.69811174 | 0.35664179 | -7.8988241 |  |
| [C5xC5](../figures/statistics/p5/sylow_C5xC5.png) | 480 | 498,293 | 0.001526734 | 0.0015840267 | 0.71823719 | 0.71737201 | 0.75612808 | -0.0027752527 |  |
| [C625](../figures/statistics/p5/sylow_C625.png) | 500 | 461,442 | 0.0014138252 | 0.0015206656 | 0.33544318 | 0.36043013 | 0.42665354 | -4.890927 |  |
| [C25xC5](../figures/statistics/p5/sylow_C25xC5.png) | 2,000 | 120,159 | 0.00036815855 | 0.0003801664 | 0.82448601 | 0.72497418 | 0.61927914 | -0.010012386 |  |
| [C3125](../figures/statistics/p5/sylow_C3125.png) | 2,500 | 52,630 | 0.00016125454 | 0.00030413312 | 0.69685816 | 0.73187342 | 0.78849947 | -0.0035304391 |  |
| [C125xC5](../figures/statistics/p5/sylow_C125xC5.png) | 10,000 | 22,487 | 6.8898554e-05 | 7.603328e-05 | 0.48238528 | 0.47948062 | 0.51336573 | -0.053820018 |  |
| [C15625](../figures/statistics/p5/sylow_C15625.png) | 12,500 | 2,047 | 6.2718611e-06 | 6.0826624e-05 | 0.92832061 | 0.94573981 | 0.96075055 | -3.7506972e-05 | Sparse |
| [C625xC5](../figures/statistics/p5/sylow_C625xC5.png) | 50,000 | 2,562 | 7.8497841e-06 | 1.5206656e-05 | 0.68495275 | 0.76039393 | 0.80569348 | -0.00012713858 | Sparse |
| [C3125xC5](../figures/statistics/p5/sylow_C3125xC5.png) | 250,000 | 94 | 2.8800925e-07 | 3.0413312e-06 | 0.91884283 | 0.94589046 | 0.96411759 | -1.7649359e-06 | Sparse |
| [C25xC25](../figures/statistics/p5/sylow_C25xC25.png) | 300,000 | 744 | 2.2795626e-06 | 2.5344427e-06 | 0.45719658 | 0.47942968 | 0.37148466 | -0.036746514 | Sparse |
| [C125xC25](../figures/statistics/p5/sylow_C125xC25.png) | 1,250,000 | 100 | 3.0639282e-07 | 6.0826624e-07 | 0.73337023 | 0.73367326 | 0.81070566 | -4.7005052e-06 | Sparse |
| [C5xC5xC5](../figures/statistics/p5/sylow_C5xC5xC5.png) | 1,488,000 | 162 | 4.9635637e-07 | 5.1097634e-07 | 3.6285022 | — | 0.11978127 | -0.3953066 | Sparse, Opposite signs for gap 2 |
| [C25xC5xC5](../figures/statistics/p5/sylow_C25xC5xC5.png) | 6,000,000 | 38 | 1.1642927e-07 | 1.2672213e-07 | — | 0.66306316 | 0.27093537 | -0.012007894 | Sparse, Opposite signs for gap 1 |
| [C625xC25](../figures/statistics/p5/sylow_C625xC25.png) | 6,250,000 | 4 | 1.2255713e-08 | 1.2165325e-07 | 1 | 0.92340269 | 0.96170134 | -7.3739529e-08 | Sparse |
| [C125xC5xC5](../figures/statistics/p5/sylow_C125xC5xC5.png) | 30,000,000 | 4 | 1.2255713e-08 | 2.5344427e-08 | 1 | 1 | 0.76166367 | -5.6508002e-07 | Sparse |
| [C125xC125](../figures/statistics/p5/sylow_C125xC125.png) | 187,500,000 | 1 | 3.0639282e-09 | 4.0551082e-09 | -1.0325218 | -0.016260879 | 0.49186956 | -1.169093e-05 | Sparse |

## Sylow 7

14 types in the data. [Compare the common groups](../figures/statistics/p7/sylow_comparison.png).

| Group | Aut size | Count | Observed / $N(X)$ | CL probability | Gap 1 exponent | Gap 2 exponent | Gap 4 exponent | Gap 4 coefficient | Status |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| [trivial](../figures/statistics/p7/sylow_trivial.png) | 1 | 273,316,093 | 0.83742089 | 0.83679541 | 0.76119408 | 0.69914411 | 0.70328823 | 0.090919519 |  |
| [C7](../figures/statistics/p7/sylow_C7.png) | 6 | 45,439,736 | 0.13922409 | 0.1394659 | 1.0480283 | 0.91605235 | 0.90342122 | -0.00054761152 |  |
| [C49](../figures/statistics/p7/sylow_C49.png) | 42 | 6,489,358 | 0.019882927 | 0.0199237 | 0.85489716 | 0.66921557 | 0.82593966 | -0.00046248719 |  |
| [C343](../figures/statistics/p7/sylow_C343.png) | 294 | 899,449 | 0.0027558472 | 0.0028462429 | 0.19243134 | 0.22565241 | 0.31857043 | -39.164309 |  |
| [C7xC7](../figures/statistics/p7/sylow_C7xC7.png) | 2,016 | 132,371 | 0.00040557524 | 0.00041507709 | 0.87571743 | 0.79645728 | 0.6419763 | -0.0049420163 |  |
| [C2401](../figures/statistics/p7/sylow_C2401.png) | 2,058 | 77,139 | 0.00023634836 | 0.00040660613 | 0.68479859 | 0.71648127 | 0.76760147 | -0.0064967513 |  |
| [C49xC7](../figures/statistics/p7/sylow_C49xC7.png) | 12,348 | 20,941 | 6.4161721e-05 | 6.7767688e-05 | 0.15918962 | 0.34899986 | 0.41408162 | -0.21439411 |  |
| [C16807](../figures/statistics/p7/sylow_C16807.png) | 14,406 | 1,533 | 4.6970019e-06 | 5.8086589e-05 | 0.95220083 | 0.95923871 | 0.96958832 | -3.0543791e-05 | Sparse |
| [C343xC7](../figures/statistics/p7/sylow_C343xC7.png) | 86,436 | 1,728 | 5.294468e-06 | 9.6810982e-06 | 0.66677005 | 0.71557481 | 0.79156775 | -0.00010169118 | Sparse |
| [C2401xC7](../figures/statistics/p7/sylow_C2401xC7.png) | 605,052 | 25 | 7.6598205e-08 | 1.383014e-06 | 0.95667115 | 0.97179842 | 0.97944959 | -6.0882378e-07 | Sparse |
| [C49xC49](../figures/statistics/p7/sylow_C49xC49.png) | 4,840,416 | 27 | 8.2726062e-08 | 1.7287675e-07 | 0.40564141 | 0.6406282 | 0.88544277 | -2.9671903e-07 | Sparse |
| [C343xC49](../figures/statistics/p7/sylow_C343xC49.png) | 29,647,548 | 1 | 3.0639282e-09 | 2.8224776e-08 | 1.1873353 | 0.91710919 | 0.9585546 | -1.8106565e-08 | Sparse |
| [C7xC7xC7](../figures/statistics/p7/sylow_C7xC7xC7.png) | 33,784,128 | 2 | 6.1278564e-09 | 2.4768892e-08 | 1 | 0.79497632 | 0.89748816 | -4.7760153e-08 | Sparse |
| [C49xC7xC7](../figures/statistics/p7/sylow_C49xC7xC7.png) | 203,297,472 | 1 | 3.0639282e-09 | 4.1161132e-09 | -0.96789463 | 0.016052684 | 0.50802634 | -8.8690755e-06 | Sparse |

[All dyadic checkpoints and pairs](dyadic_estimates.csv) · [Endpoint estimates](fit_summary.csv) · [JSON](fits.json)
