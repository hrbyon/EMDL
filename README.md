# ATR-FTIR Kramers–Krönig 보정 도구

ATR-FTIR 스펙트럼을 투과형 FTIR과 비교 가능한 흡광 스펙트럼으로 보정합니다.
보정 절차는 다음 논문을 그대로 구현한 것입니다.

> Taewoo Kim, Minhaeng Cho, Kyungwon Kwak,
> *Quantitative Analysis of the Li-Ion Solvation Structure in Li-Ion Battery Electrolytes Using ATR-FTIR Spectroscopy*,
> **Anal. Chem.** 2024, 96, 15924–15930. [doi:10.1021/acs.analchem.4c02779](https://doi.org/10.1021/acs.analchem.4c02779)

## 왜 필요한가

ATR-FTIR은 반사를 측정하므로 Beer–Lambert 법칙이 그대로 성립하지 않습니다.
흡광이 강한 시료(예: 전해액)에서는 Kramers–Krönig 관계에 따라 흡수 피크 근처에서 굴절률이 크게 요동하고,
그 결과 증발파의 침투 깊이가 파수에 따라 달라져 **피크가 red shift 하고 세기가 왜곡**됩니다.
장비에 기본 탑재된 ATR 보정(침투 깊이 정규화)은 흡광이 약한 시료에만 유효합니다.

## 보정 절차

1. pATR → 반사율 `R` → `Rs` (45°에서 Abelès 관계 `Rp = Rs²`)
2. `ln√Rs`에 대한 수정 Kramers–Krönig 적분(Maclaurin 방식, Ohta & Ishida 1988)으로 위상 `Θs` 산출
3. Fresnel 식 역산으로 복소굴절률 `n̂ = n + ik` 계산
4. `A = 1.74πkνl` 로 흡광 스펙트럼 재구성 (논문 식 16)

`Θs(∞)`는 시료의 고주파 굴절률 `n(∞)`로 결정되므로, 이 값을 별도로 정하거나 추정해야 합니다.

## 구성

| 파일 | 내용 |
| --- | --- |
| `atr_kk.py` | 보정 라이브러리. 단독 실행하면 모사 스펙트럼으로 자체 검증을 수행합니다. |
| `cof_xrd.py` | COF 분말 XRD 분석 라이브러리·명령줄 도구. `python cof_xrd.py --selftest`로 검증합니다. |
| `atr-workbench.html` | 브라우저에서 파일을 올려 보정하고 그래프·CSV를 얻는 단일 파일 웹 도구. 외부 의존성이 없습니다. |

## 사용법

### 1. Python

```python
import numpy as np
from atr_kk import atr_kk_correct

nu = ...   # 파수 (cm-1), 오름차순, 등간격
R  = ...   # 측정 반사율 (0-1, 비편광)

out = atr_kk_correct(nu, R, n_inf=1.40, n0=2.4, theta_deg=45.0, path_um=1.0)
out["A_corr"], out["n"], out["k"]
```

자체 검증 실행:

```bash
python atr_kk.py
```

알려진 Lorentz 진동자로 ATR 스펙트럼을 모사한 뒤 역산합니다.
원시 ATR에서 6–7 cm⁻¹ red shift 된 밴드를 **0.1 cm⁻¹ 이내**로 복원하며, `k` 오차는 6% 이내입니다.

### 2. 웹 워크벤치

`atr-workbench.html`을 브라우저로 열면 됩니다. 서버가 필요 없습니다.

- 입력: `.csv .txt .dat .tsv` 또는 붙여넣기. 쉼표·탭·공백 구분과, `x + 다중 y` / `(x,y)(x,y)…` 쌍 배열을 자동 인식합니다.
- y축: %T/%R, T/R, 흡광도, 직접 지정한 전스케일 값 중 선택. 기준선 구간의 중앙값을 100%로 삼습니다.
- ATR 조건: 결정 n₀(ZnSe·다이아몬드·Ge·KRS-5·직접 입력), 입사각, n(∞)(시료별 개별 지정 가능), 경로 길이
- 출력: 원시/보정 겹쳐 보기, `n(ν)`, `k(ν)`, 부분보정 비교, 구간별 피크 표(원시·보정·Δ), CSV, 그래프 PNG
- `n(∞)` 스캔으로 배위수·피크 위치의 민감도를 확인할 수 있습니다.

파수 간격이 일정하지 않으면 균일 격자로 보간합니다. KK 적분에 균일 격자가 필요하기 때문입니다.

## 결과를 읽을 때 주의할 점

- **R_min이 0.02 아래**로 내려간 밴드는 전반사 조건이 깨진 것입니다. 위치도 세기도 신뢰할 수 없습니다. 워크벤치가 경고를 띄웁니다.
- **n(∞)를 잘못 잡으면** 배위수가 과대·과소평가됩니다. 논문에서도 1 M DMC/LiClO₄의 CN이 n(∞)에 따라 3.05–3.2로 흔들렸습니다.
- **부분보정**(장비 기본 ATR correction)은 흡광이 강한 전해액에서 거의 효과가 없습니다. 논문 기준 CN 4.0 → 3.9.
- KK 적분 구간이 제한되고 이상적인 유전체 계면을 가정하므로, 보정 후에도 투과 FTIR과 완전히 일치하지는 않습니다.

## 인용

이 코드를 사용해 얻은 결과를 발표할 때는 위 논문을 인용하십시오.

---

# COF PXRD 워크벤치 (`xrd-workbench.html`)

COF 분말 XRD의 피크 분석, Pawley 정제, 시뮬레이션·구조 모델 비교를 브라우저에서 수행합니다. 외부 의존성이 없는 단일 HTML 파일입니다.

페이지는 세 부분입니다. **A 구조 시뮬레이션**은 구조 모델의 층을 AA·inclined·serrated·AB로 쌓고 치환기 정렬 비율을 바꿔 가며 원자 고정 fit으로 실험과 비교합니다. **B Pawley 정제**는 원자 좌표 없이 격자·강도·결정 크기를 구합니다. **C 종합 결론**은 두 결과를 합쳐 격자·적층·층간거리·치환기 배치를 판정합니다.

## 구조 시뮬레이션 (A)

- **라이브러리**: xtb 6.7.1 GFN-FF(주기 경계)로 최적화한 층 — TpBd-(SO3Li)2, TpBD, TpPa-SO3Li, 2D-PAI, COF-1. CIF를 올리면 첫 층을 잘라 씁니다.
- **적층 스캔**: AA, AB, inclined·serrated(1–4 Å, 0°/30°)를 모두 fit해 Rwp로 순위를 매기고, 계산해 둔 GFN-FF 적층 에너지를 옆에 표시합니다.
- **무질서 스캔**: 지정한 치환기 원소(예: `S,Li`, 술포네이트 O 포함)의 점유율을 0–1로 바꿔 fit합니다.
- **내보내기**: 현재 적층 모델의 CIF, xtb/DFT용 POSCAR.
- **신뢰도**: 문헌 DFT/DFTB 구조로 검증한 결과 면내 격자 ±2 %, 층간거리 0.15–0.35 Å 과소, 중성 골격의 적층 에너지 순서는 문헌과 일치합니다. Li+를 포함한 계는 GFN-FF 에너지가 불연속이라 적층 에너지를 쓰지 않습니다. 벤치마크 표는 페이지의 "계산 신뢰도" 항목에 있습니다.

## Pawley 정제 (B)

- **입력**: `(2θ, I)(2θ, I)…` 쌍 배열(시뮬레이션·실험 열 길이가 달라도 됨) 또는 `2θ + 다중 I`. 시뮬레이션 열은 자동 인식합니다. 선택적으로 CIF(대칭 연산자·점유율 지원)를 올려 모델 강도를 계산합니다.
- **피크 분석**: 피크 위치·d·FWHM·Scherrer 결정 크기, 육방정/정방정 자동 지수화, π–π 적층 피크와 층 수 추정.
- **Pawley 정제**: 격자상수(a, c), zero shift, Thompson–Cox–Hastings pseudo-Voigt 폭(면내·적층 방향 분리), 비대칭, 저각 지수 배경 + Chebyshev 배경. 비선형 파라미터는 Levenberg–Marquardt, 반사 강도와 배경은 NNLS로 풉니다.
- **LP 가중 프로파일**: 5° 이하의 넓은 피크는 Lorentz–편광 인자가 피크 폭 안에서 크게 변해 저각으로 치우칩니다. 이를 프로파일에 반영합니다.
- **zero 스캔**: zero shift를 ±0.15°로 고정하며 격자상수·Rwp 변화를 보여 줍니다. 저각 반사만 있는 COF에서 a와 zero의 상관을 확인하는 데 필요합니다.
- **모델 비교**: Pawley로 추출한 상대 강도를 시뮬레이션·CIF 강도와 비교해 3σ 이상 벗어나는 반사를 표시하고, 모델 강도를 고정한 fit의 Rwp를 계산합니다.
- **자동 해석**: 셀 불일치, 결정 크기, 잔차 상관(Durbin–Watson), 강도 불일치, a–zero 상관을 요약합니다.
- **출력**: 패턴 CSV(obs/calc/bg/diff), 반사 표 CSV, 그래프 PNG, 요약 텍스트.

## Python (`cof_xrd.py`)

웹 워크벤치와 같은 절차를 스크립트로 수행합니다. 필요한 패키지는 `numpy`, `scipy`이고, 그림을 저장할 때만 `matplotlib`을 씁니다.

```bash
python cof_xrd.py data.csv                              # 실험·시뮬레이션 열 자동 인식, Pawley, 시뮬레이션 비교
python cof_xrd.py data.csv --cif model.cif --scan --plot --out result
python cof_xrd.py --selftest                            # 합성 육방정 COF로 a·강도 복원 검증
```

주요 옵션: `--lattice hex|tet`, `--range LO HI`(정제 구간), `--zero 0.0`(zero 고정), `--exp/--sim`(열 번호 지정), `--lam`(파장), `--B`(CIF 강도용 등방 온도인자).
출력: `*_pawley_fit.csv`(obs/calc/bg/diff), `*_reflections.csv`(반사별 강도·±σ·모델 강도), `*_pawley.png`.

```python
import cof_xrd as cx
series = cx.load_xrd("data.csv")
exp = next(s for s in series if not s.is_sim)
fit = cx.pawley(exp.x, exp.y, lattice="hex", tt_range=(2.3, 12))
print(fit.summary())
scan = cx.zero_scan(exp.x, exp.y, lattice="hex", tt_range=(2.3, 12))
st = cx.read_cif("model.cif")
_, d, tt, I = cx.structure_reflections(st, tt_max=20)
pct, _ = cx.model_percentages(fit, d, I)
print(cx.compare(fit, pct))                             # 3σ 이상 다른 반사
```

## 불확도

±σ는 공분산 행렬에서 구하고 잔차의 직렬 상관에 따라 √(2/DW)배 키웁니다. zero shift 상관에 의한 계통 오차는 포함되지 않으므로 zero 스캔 결과를 함께 보고하세요.

## 참고 문헌

- G. S. Pawley, *J. Appl. Cryst.* **1981**, 14, 357–361.
- P. Thompson, D. E. Cox, J. B. Hastings, *J. Appl. Cryst.* **1987**, 20, 79–83.
- D. T. Cromer, J. B. Mann, *Acta Cryst.* **1968**, A24, 321–324 (원자 산란인자).
