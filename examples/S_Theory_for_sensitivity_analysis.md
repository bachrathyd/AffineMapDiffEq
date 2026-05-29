A módszer pontos matematikai leírása
Jelölések
$\mathcal{M}(p) : X \to X$ — a monodrómia-operátor (periódus-leképezés), $X$ a késő-argumentumú rendszer véges-dimenziós diszkretizált históriájának tere
$p = (p_1, \ldots, p_n) \in \mathbb{R}^n$ — paramétervektor
$\mu_j(p)$ — Floquet-multiplikátorok (sajátértékek), $|\mu_j| < 1$ ⟺ stabil
1. lépés — Nominális Krylov–Schur bázis (Float64)
A Krylov-módszer (jelen esetben KrylovKit.schursolve) megadja a domináns $n_b$ sajátpárt a nominális $p_0$-nál:

$$\mathcal{M}(p_0), q_j \approx \mu_j^{(0)}, q_j, \qquad j = 1, \ldots, n_b$$

ahol ${q_j}$ közel-ortonormált Schur-bázis. Ez Float64 aritmetikával fut — ide nem kerül Dual szám.

2. lépés — Dual-számos paraméter-seed
Az $n_{active}$ aktív (nem nulla bizonytalanságú) paramétert egyszerre seedeljük, chunked forward-mode AD-dal:

$$\tilde{p}i = p_i^{(0)} + \delta{i,\text{active}(k)},\varepsilon_k, \qquad k = 1, \ldots, n_{active}$$

Julia-típus: Dual{SensTag, Float64, n_active}, azaz minden szám egy értéket és $n_{active}$ darab parciálist hordoz. Egyszerre seedelünk, így az összes $\partial/\partial p_i$ egyszerre "áramlik" végig az integrátoron.

3. lépés — Tömörített operátor felépítése
A $Q = [q_1 | \cdots | q_{n_b}]$ bázisra vetítve az operátort $n_b$ LinMap-hívással:

$$\widetilde{H}_{ab} = \langle q_a,, \mathcal{M}(\tilde{p}), q_b \rangle, \qquad a, b = 1, \ldots, n_b$$

Minden $\widetilde{H}_{ab}$ egy Dual{SensTag} szám:

$$\widetilde{H}{ab} = \underbrace{H{ab}(p_0)}{\text{Float64 érték}} + \sum{i=1}^{n_{active}} \underbrace{\frac{\partial H_{ab}}{\partial p_i}}_{\text{partial}} \varepsilon_i$$

ahol:
$$\frac{\partial H_{ab}}{\partial p_i} = \left\langle q_a,, \frac{\partial \mathcal{M}}{\partial p_i} q_b \right\rangle$$

A Dual-aritmetika automatikusan kiszámolja ezeket — nem kell expliciten deriválni az integrátort.

A LinMap-hívás belsejében a belső AffineTag-es Dual réteg (ami a $q_b$ irányba perturbál) fölé kerül a SensTag-es Dual réteg (ami a paraméterekben perturbál). A partialpart az AffineTag réteg deriváltját hántja le, a SensTag partialok bennmaradnak az eredményben.

4. lépés — Első rendű sajáréték-perturbáció az $n_b \times n_b$ mátrixon
Az értékrész adja a nominális kis mátrixot:

$$H_0 = \text{value}!\left(\widetilde{H}\right) \in \mathbb{R}^{n_b \times n_b}$$

Ezt Float64 eigen-nel felbontjuk:

$$H_0, R = R, \text{diag}(\lambda_1, \ldots, \lambda_{n_b}), \qquad L^\top = R^{-1}$$

ahol $r_m, \ell_m$ az $m$-edik jobb- ill. bal-sajátvektor. Ezután a $p_i$ szerinti perturbációmátrix:

$$H_1^{(i)} = \text{partials}_i!\left(\widetilde{H}\right) \in \mathbb{R}^{n_b \times n_b}$$

Az első rendű sajáréték-perturbáció (Hellmann–Feynman-tétel véges dimenziós formája):

$$\frac{\partial \lambda_m}{\partial p_i} = \left(L^\top H_1^{(i)}, R\right)_{mm} = \ell_m^\top \frac{\partial H}{\partial p_i} r_m$$

5. lépés — $|\mu_m|$ érzékenysége
A lánc-szabállyal:

$$\frac{\partial |\mu_m|}{\partial p_i} = \text{Re}!\left(\frac{\overline{\mu_m}}{|\mu_m|} \cdot \frac{\partial \mu_m}{\partial p_i}\right) = \frac{\text{Re}!\left(\overline{\lambda_m} \cdot (L^\top H_1^{(i)} R)_{mm}\right)}{|\lambda_m|}$$

Ez az implementáció (robust_control.jl:237-241):


dabs[m, i] = real(conj(λ[m]) * M1[m, m]) / am   # M1 = Lt * H1 * R
Összes lépés összefoglalva
$$\boxed{\frac{\partial |\mu_m|}{\partial p_i} \approx \frac{\text{Re}!\left(\overline{\lambda_m} \cdot \ell_m^\top \left\langle Q^\top, \frac{\partial \mathcal{M}}{\partial p_i} Q\right\rangle_{:,m}\right)}{|\lambda_m|}}$$

ahol $Q = [q_1|\cdots|q_{n_b}]$, $\lambda_m$ és $\ell_m, r_m$ az $n_b \times n_b$ tömörített $H_0$ sajátpárjai.

Közelítés forrásai (és hogyan csökken $n_b$ növelésével)
1. A Schur-bázis csonkítása: A bal-sajátvektor $\ell_m$ (ami az érzékenységben megjelenik) csak akkor kapható pontosan, ha $\ell_m \in \text{span}(Q)$. Ha $n_b$ kicsi, a teljes bal-sajátvektort $Q$ nem fedi le → az érzékenység torzított. Növelve $n_b$-t, a Krylov-altér kibővül és $\ell_m$ egyre pontosabban benne van.

2. Első rendű perturbáció: Az összefüggés $\lambda_m(p) \approx \lambda_m(p_0) + \nabla_p \lambda_m \cdot \Delta p$ csak kis $\Delta p$-re egzakt. A módszer elvileg pontosan adja a loká­lis deriváltat — nincsen FD lépéshiba.

Miért O(n) és nem O(n²)?
A szokásos forward-mode AD érve: az $n_{active}$ partial egy integrációs menetben áramlik végig az integrátoron. Egy összeadás/szorzás $O(n_{active})$ FLOP-ba kerül, de nem kell $n_{active}$-szor újraindítani az integrációt. Ezért a teljes költség:

$$\underbrace{n_b}{\text{LinMap hívások}} \times \underbrace{O(n{active})}{\text{Dual aritmetika/lépés}} \times T{\text{integrate}}$$

Mivel $n_b$ rögzített algoritmikus paraméter: O(n) a paraméterszámban. A $2^n$ sarokpontos brute-force bejárás ezzel szemben exponenciális — 5 paraméternél már 32, 10-nél 1024 teljes MDBM futást igényelne