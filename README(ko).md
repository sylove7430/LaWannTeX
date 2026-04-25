# LaWannTeX

## 실행 방법 및 필수 파일들
LaWannTeX는 Wannier90 Hamiltonian를 읽고 추출한 호핑 파라미터를 TeX 파일로 정리해주는 Julia 프로젝트입니다.
루트에 있는 `main.jl`이 이 프로그램의 엔트리포인트이고, 내부 구현은 `src/` 아래에 정리되어 있습니다. 실행하는 방법으로는 아래와 같습니다: 

```bash
julia main.jl --input examples/graphene/input.toml
```

이때 입력 파일은 top-level 태그인 `[files]`를 포함한 TOML (input.toml)만 있으면 됩니다. 필요한 파일들에는 Wannier90 input file인 `wannier90.win`와 Wannier90 output hr file인 `wannier90_hr.dat`이 있습니다. `.win` 파일에서는 lattice와 projection 라벨을 읽습니다. 손으로 직접 쓰려면 같은 TOML 안에 top-level `[structure]`와 `[[orbitals]]`를 넣으면 됩니다. `_hr.dat` 파일에서는 tight-binding parameter를 읽습니다.

## 입력 파일 `.toml` 작성 방법

- `[files]`에는 읽어들일 Wannier90 Hamiltonian 파일 `hr`과 tex 파일을 저장할 곳인 `out`이 필요합니다.
- `[render]`에는 보통 `atol`, `include_second_quant`를 넣습니다. 특정 오비탈을 선택하지 않고 제일 큰 호핑을 정렬 후에 뽑는다면 `max_rblocks`을 줄 수 있습니다. 이들은 각각 호핑 파라미터의 허용 오차, second quantization으로 나타낼 것인지, 그리고 $\mathbf{R}$를 어디까지 뽑을 것인지를 결정합니다.
- `[selection]`이 없으면 기본적으로 정렬 후 작성 모드(legacy)로 동작합니다.
- `[selection]`에서 `mode = "orbital"`을 주면 선택한 Wannier index 집합끼리의 호핑 파라미터만 읽고, `max_nn` shell까지의 항만 남깁니다.
- 이때 `max_nn` shell은 [Kwant](https://kwant2.uber.space)의 `neighbors(n)`처럼 selected basis 전체에서 global non-zero distance class로 계산합니다. `max_nn = 0`은 raw home-cell zero-distance onsite 항만 남깁니다.
- 구조와 오비탈 정보는 `wannier90.win`에서 읽거나, `[structure]`와 `[[orbitals]]`로 직접 적는 두 방식 중 하나를 쓰면 됩니다.
- `[[orbitals]].orbitals`와 `wannier90.win` projection에서는 shorthand를 지원합니다. 현재 `s`, `p`, `d`, `f`가 Wannier90 순서의 real orbital label로 자동 확장됩니다.

## 예제

예제 파일은 examples 폴더에서 확인할 수 있습니다.

### 1. `wannier90.win`을 사용하는 경우:

```toml
[files]
hr = "wannier90_hr.dat"
win = "wannier90.win"
out = "reports/hamiltonian.tex"

[render]
atol = 1.0e-6
max_rblocks = 12
include_second_quant = true

[selection]
mode = "legacy"
```

### 2. 구조와 오비탈을 직접 적는 경우:

```toml
[structure]
lattice = [
  [2.46, 0.0, 0.0],
  [1.23, 2.130422493, 0.0],
  [0.0, 0.0, 20.0],
]

[[orbitals]]
name = "C"
position_frac = [0.0, 0.0, 0.0]
orbitals = ["pz"]
spins = ["none"]

[[orbitals]]
name = "C"
position_frac = [0.3333333333, 0.3333333333, 0.0]
orbitals = ["pz"]
spins = ["none"]

[files]
hr = "wannier90_hr.dat"
out = "reports/hamiltonian.tex"

[render]
atol = 1.0e-6
max_rblocks = 12
include_second_quant = true

[selection]
mode = "legacy"
```

### 3. 오비탈 선택 모드(`orbital`)를 사용하는 경우:

```toml
[files]
hr = "wannier90_hr.dat"
out = "reports/orbital_model.tex"

[render]
atol = 1.0e-6
include_second_quant = true

[selection]
mode = "orbital"
wannier_indices = [1, 2]
max_nn = 1

[structure]
lattice = [
  [2.46, 0.0, 0.0],
  [1.23, 2.130422493, 0.0],
  [0.0, 0.0, 20.0],
]

[[orbitals]]
name = "C"
position_frac = [0.0, 0.0, 0.0]
orbitals = ["pz"]
spins = ["none"]

[[orbitals]]
name = "C"
position_frac = [0.3333333333, 0.3333333333, 0.0]
orbitals = ["pz"]
spins = ["none"]
```
