#!/usr/bin/env python3
"""Fix aggregate mean_slope/mean_r2 in calibration_validation.json."""
import json
from pathlib import Path

cal_path = Path(__file__).parent.parent.parent / "benchmark" / "calibration_validation.json"
cal = json.loads(cal_path.read_text())
vals = cal["validations"]

mean_slope = sum(v["calibration_slope"] for v in vals) / len(vals)
mean_r2 = sum(v["r_squared"] for v in vals) / len(vals)
well = sum(1 for v in vals if v["is_well_calibrated"])

print(f"mean_slope: {cal['mean_slope']} -> {mean_slope:.4f}")
print(f"mean_r2: {cal['mean_r2']} -> {mean_r2:.4f}")
print(f"well: {cal['well_calibrated']} -> {well}/{len(vals)}")

cal["mean_slope"] = round(mean_slope, 4)
cal["mean_r2"] = round(mean_r2, 4)
cal["well_calibrated"] = well
cal["poorly_calibrated"] = len(vals) - well
cal["global_status"] = "GOOD" if well > len(vals) * 0.7 else "NEEDS_IMPROVEMENT"

cal_path.write_text(json.dumps(cal, indent=2))
score = min(100, max(10, mean_r2 * 100 - abs(mean_slope - 1) * 50))
print(f"Calibration score: {score:.1f}/100")
