#!/usr/bin/env bash
# Usage: bash run_wsl.sh OUT_DIRECTORY [smoke|grid|paper]
set -euo pipefail

source_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: bash run_wsl.sh OUT_DIRECTORY [smoke|grid|paper]" >&2
    exit 2
fi
mode=${2:-smoke}
case "$mode" in
    smoke|grid|paper) ;;
    *) echo "Expected smoke, grid or paper" >&2; exit 2 ;;
esac
launcher=${VASP_LAUNCHER:-vasp-cpu}
potcars=${VASP_POTCAR_ROOT:?Set VASP_POTCAR_ROOT to your local licensed PBE PAW root}
mn=${VASP_MN_PAW:-Mn_pv}
ranks=${VASP_RANKS:-4}
limit=${VASP_TIMEOUT_SECONDS:-480}
if ! executable=$(command -v -- "$launcher"); then
    echo "Set VASP_LAUNCHER to a launcher accepting MPI ranks and VASP variant" >&2
    exit 2
fi
launcher=$executable
if [[ ! -r "$potcars/$mn/POTCAR" || ! -r "$potcars/F/POTCAR" ]]; then
    echo "Expected readable PAW files for $mn and F under VASP_POTCAR_ROOT" >&2
    exit 2
fi
mkdir -- "$1"
output=$(cd -- "$1" && pwd)
scratch=$(mktemp -d /tmp/ammoire-mnf2.XXXXXXXX)
trap 'rm -f -- "$scratch/scf/POTCAR" "$scratch/probe/POTCAR"' EXIT
mkdir -- "$scratch/scf" "$scratch/probe"

printf 'mode=%s\nmpi_ranks=%s\nPAW=%s,F\nscratch=%s\n' "$mode" "$ranks" "$mn" "$scratch" > "$output/run.txt"
{
    sha256sum "$potcars/$mn/POTCAR" "$potcars/F/POTCAR" | awk '{print $1}'
    grep -m1 'TITEL' "$potcars/$mn/POTCAR"
    grep -m1 'TITEL' "$potcars/F/POTCAR"
} > "$output/paw.txt"
cat -- "$potcars/$mn/POTCAR" "$potcars/F/POTCAR" > "$scratch/scf/POTCAR"
ln -s -- "$scratch/scf/POTCAR" "$scratch/probe/POTCAR"
cp -- "$source_dir/POSCAR" "$scratch/scf/POSCAR"
cp -- "$source_dir/POSCAR" "$scratch/probe/POSCAR"
cp -- "$source_dir/INCAR.scf" "$scratch/scf/INCAR"
cp -- "$source_dir/INCAR.probe" "$scratch/probe/INCAR"
cp -- "$source_dir/KPOINTS.$mode" "$scratch/scf/KPOINTS"
cp -- "$source_dir/KPOINTS.probe" "$scratch/probe/KPOINTS"

run_stage() {
    local stage=$1 code=0 file
    mkdir -- "$output/$stage"
    echo "Starting VASP $stage ($mode, $ranks MPI ranks)"
    (cd -- "$scratch/$stage"; timeout --kill-after=20 "$limit" "$launcher" "$ranks" std) \
        > "$output/$stage/vasp.log" 2>&1 || code=$?
    for file in INCAR POSCAR KPOINTS OUTCAR OSZICAR EIGENVAL PROCAR; do
        if [[ -f "$scratch/$stage/$file" ]]; then
            cp -- "$scratch/$stage/$file" "$output/$stage/$file"
        fi
    done
    echo "$code" > "$output/$stage/exit_code.txt"
    return "$code"
}

run_stage scf
if ! grep -q 'aborting loop because EDIFF is reached' "$scratch/scf/OUTCAR" ||
   [[ ! -s "$scratch/scf/CHGCAR" ]]; then
    echo "SCF did not produce converged CHGCAR; see $output/scf" >&2
    exit 1
fi
ln -s -- "$scratch/scf/CHGCAR" "$scratch/probe/CHGCAR"
run_stage probe
if ! grep -q 'aborting loop because EDIFF is reached' "$scratch/probe/OUTCAR" ||
   [[ ! -s "$scratch/probe/EIGENVAL" ]]; then
    echo "Probe did not produce converged eigenvalues; see $output/probe" >&2
    exit 1
fi
echo "Finished $output"
