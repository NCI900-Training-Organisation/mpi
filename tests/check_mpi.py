#!/usr/bin/env python3
"""Build and check the workshop reference solutions using a local MPI runtime."""
import math
import os
from pathlib import Path
import re
import shutil
import signal
import struct
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = {
    "day1": ("blocking", "nonblocking", "persistent"),
    "day2": ("collective", "win", "io"),
}


def execute(command, cwd, success=True):
    # A separate process group lets a timeout clean up the launcher and its ranks.
    process = subprocess.Popen(command, cwd=cwd, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True,
                               start_new_session=True)
    try:
        output, _ = process.communicate(timeout=30)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        output, _ = process.communicate()
        raise AssertionError(f"Timed out: {command}\n{output}") from None
    assert (process.returncode == 0) == success, (command, process.returncode, output)
    return output


def expected_residual(size, iterations):
    return math.pi ** 2 / (size - 1) * math.cos(math.pi / (size - 1)) ** iterations


def check_output(path, size, iterations):
    # The sine forcing is a discrete eigenvector, so every Jacobi iterate is known.
    data = path.read_bytes()
    assert len(data) == size * size * struct.calcsize("d"), (size, len(data))
    values = struct.unpack(f"{size * size}d", data)
    h = 1.0 / (size - 1)
    rho = math.cos(math.pi * h)
    factor = (math.pi ** 2 * h ** 2 / 2) * (1 - rho ** iterations) / (1 - rho)
    for i in range(size):
        for j in range(size):
            expected = factor * math.sin(math.pi * i * h) * math.sin(math.pi * j * h)
            actual = values[i * size + j]
            assert math.isfinite(actual) and abs(actual - expected) < 1e-11, (
                size, iterations, i, j, actual, expected)


def main():
    for tool in ("make", "mpicc", "mpiexec"):
        if shutil.which(tool) is None:
            raise SystemExit(f"Required command not found: {tool}")
    runs = 0
    with tempfile.TemporaryDirectory(prefix="mpi-regression-") as directory:
        work = Path(directory)
        binaries = {}
        for day, variants in VARIANTS.items():
            build = work / day
            shutil.copytree(ROOT / day / "solution", build)
            # Use the actual Makefiles; source-only copies prevent stale executables.
            for filename in build.iterdir():
                if filename.suffix not in (".c", ".h") and filename.name != "Makefile":
                    if filename.is_file():
                        filename.unlink()
            execute(["make", *variants, "CFLAGS=-g -Wall -Wextra -Werror -O3"], build)
            for variant in variants:
                binaries[variant] = build / f"laplace_mpi_{variant}"

        for variant, executable in binaries.items():
            folder = work / f"run-{variant}"
            folder.mkdir()
            # One rank, multiple ranks, nondivisible interior rows, and zero iterations.
            for ranks, size, iterations in ((1, 18, 10), (2, 18, 10),
                                            (4, 19, 10), (2, 18, 0)):
                output = execute(["mpiexec", "-np", str(ranks), str(executable),
                                  str(size), str(iterations), "Jacobi"], folder)
                match = re.search(r"Final [Rr]esidual\s+(\S+)", output)
                assert match, output
                residual = float(match[1])
                assert math.isfinite(residual) and abs(
                    residual - expected_residual(size, iterations)) < 1e-6, output
                assert "leaked" not in output.lower(), output
                if variant == "io":
                    check_output(folder / "laplace-soln-whole", size, iterations)
                runs += 1
            # Missing arguments must fail promptly on every variant.
            execute(["mpiexec", "-np", "2", str(executable)], folder, success=False)
            runs += 1
            print(f"PASS {variant}: numerical results and collective argument failure", flush=True)

        for args in (("18", "-1", "Jacobi"), ("18", "10", "invalid"),
                     ("5", "10", "Jacobi"), ("18.5", "10", "Jacobi"),
                     ("invalid", "10", "Jacobi"),
                     ("999999999999999999999", "10", "Jacobi")):
            execute(["mpiexec", "-np", "2", str(binaries["blocking"]), *args],
                    work, success=False)
            runs += 1

        # Reuse an existing output file, then reduce the grid size.
        for size in (18, 10):
            execute(["mpiexec", "-np", "2", str(binaries["io"]),
                     str(size), "10", "Jacobi"], work)
            check_output(work / "laplace-soln-whole", size, 10)
            runs += 1
        print(f"PASS {runs} MPI runs, including output contents and file resizing")


if __name__ == "__main__":
    main()
