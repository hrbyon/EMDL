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
