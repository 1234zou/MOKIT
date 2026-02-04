#!/usr/bin/env bash
set -euo pipefail

bin_dir="${1:-bin}"

pass_count=0
xfail_count=0
fail_count=0

found_any=0

for test_bin in "${bin_dir}"/*; do
  if [[ ! -x "${test_bin}" ]]; then
    continue
  fi

  found_any=1

  base_name="$(basename "${test_bin}")"
  echo "=== START ${base_name} ==="

  if [[ "${base_name}" == *_xfail ]]; then
    set +e
    output="$(${test_bin} 2>&1)"
    status=$?
    set -e
    if [[ ${status} -eq 2 ]]; then
      echo "${output}"
      echo "FAIL: ${base_name} expected failure but passed"
      echo "=== END ${base_name} ==="
      fail_count=$((fail_count + 1))
    elif [[ ${status} -ne 0 ]]; then
      echo "${output}"
      echo "XFAIL: ${base_name}"
      echo "=== END ${base_name} ==="
      xfail_count=$((xfail_count + 1))
    elif echo "${output}" | grep -q "ERROR in subroutine"; then
      echo "${output}"
      echo "XFAIL: ${base_name}"
      echo "=== END ${base_name} ==="
      xfail_count=$((xfail_count + 1))
    else
      echo "${output}"
      echo "FAIL: ${base_name} expected failure but passed"
      echo "=== END ${base_name} ==="
      fail_count=$((fail_count + 1))
    fi
  else
    set +e
    output="$(${test_bin} 2>&1)"
    status=$?
    set -e
    if [[ ${status} -eq 0 ]]; then
      echo "${output}"
      echo "PASS: ${base_name}"
      echo "=== END ${base_name} ==="
      pass_count=$((pass_count + 1))
    else
      echo "${output}"
      echo "FAIL: ${base_name}"
      echo "=== END ${base_name} ==="
      fail_count=$((fail_count + 1))
    fi
  fi
done

if [[ ${found_any} -eq 0 ]]; then
  echo "No test binaries found in ${bin_dir}" 1>&2
  exit 1
fi

total_count=$((pass_count + xfail_count + fail_count))
echo "===== TEST SUMMARY ====="
echo "Summary: pass=${pass_count} xfail=${xfail_count} fail=${fail_count} total=${total_count}"
echo "========================"

if [[ ${fail_count} -ne 0 ]]; then
  echo "${fail_count} test(s) failed"
  exit 1
fi
