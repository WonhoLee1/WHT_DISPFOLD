# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Project: WHT_DISPFOLD — 폴더블 디스플레이 폴딩 해석 Pre-Processor

FEM 멀티-솔버(Multi-Solver) 폴더블 디스플레이 적층 구조 폴딩 해석 전처리기.

5층 A-B-A-B-A 적층 구조, 힌지(35-65mm) 기준 90° 폴딩을 **4개 솔버**에서 동일한 메시/물성으로 해석할 수 있는 입력 파일 생성.

### 지원 솔버

| 솔버 | 요소 | 파일 | 설명 |
|------|------|------|------|
| **OptiStruct 2D** | CQPSTN + PPLANE (plane strain) | `display_folding_2d.fem` | 2D, 가장 가벼움 |
| **OptiStruct 3D** | CHEXA + PSOLID (solid) | `display_folding_3d.fem` | 3D solid 검증 |
| **CalculiX 3D** | C3D8R + *HYPERELASTIC | `display_folding.inp` | 오픈소스 검증 |
| **FEBio 3D** | hex8 + neo-Hookean/viscoelastic | `display_folding.feb` | 점탄성 포함 |

### 적층 구조

```
A (Neo-Hookean, E=2000MPa, nu=0.34)  — 구조 지지층, 30μm
B (Neo-Hookean, E=1~50MPa, nu=0.45) — PSA 점탄성층, 30μm
A / B / A / B / A  (5층 교번 적층)
```

### 아키텍처 의존 방향 (단방향)

```
src/display_folding.py (CLI entry)
  → src/mesh.py          (GMSH + Python fallback)
  → src/materials.py     (Neo-Hookean constants)
  → src/kinematics.py    (rotation displacement)
  → src/writers/*.py     (solver-specific output)
```

### Python 실행 환경

conda 환경: **`vdmc`**

```bash
# 올바른 실행 방법
"C:/Users/GOODMAN/miniconda3/envs/vdmc/python.exe" -m src --solver ccx

# GMSH 강제 사용
"C:/Users/GOODMAN/miniconda3/envs/vdmc/python.exe" -m src --solver opti2d --gmsh
```

### 주요 실행 명령

```bash
# CalculiX 3D 입력 생성
python -m src --solver ccx

# OptiStruct 2D (plane strain) 입력 생성
python -m src --solver opti2d

# OptiStruct 3D 입력 생성
python -m src --solver opti3d

# FEBio 입력 생성
python -m src --solver febio

# 전부 한 번에 생성
python scripts/run_all.py
```

### 개발 규칙

- **언어**: 모든 대화와 응답은 **한국어**
- **인코딩**: 모든 `open()` 호출에 `encoding='utf-8'` 명시
- **새 파일**: BOM 없는 UTF-8로 저장
- **의존성**: numpy 필수, GMSH는 옵션 (fallback 내장)
- **작업 로그**: `dev_log/issue_tracker.md` 참조

### OMO 철학

1. **단방향 의존성** — 하위 모듈이 상위 모듈을 import하지 않음
2. **외과적 수정** — 요청과 무관한 코드 개선 금지
3. **검증 가능한 완료 기준** — "작동하게 만들기"보다 구체적 조건 명시
4. **결과물 중심** — 실행 결과(.fem, .inp, .feb)가 실제 솔버에서 동작해야 함
