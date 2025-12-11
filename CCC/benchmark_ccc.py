#!/usr/bin/env python3
"""
Quick benchmark of CCC on dev set proteins.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from CCC import CCCPredictor

# Dev set proteins with expected domains
DEV_SET = {
    "A0A2T1DI41": ["62-181", "225-323", "359-452", "480-581", "678-832", "846-993", "1014-1082"],
    "A9WGU0": ["5-211", "234-331", "359-464", "490-494", "561-711", "730-867", "887-946"],
    "A0A482W0X4": ["2-324"],
    "A9WIN4": ["19-179_273-325", "217-253", "329-420"],  # discontinuous
    "P12269": ["20-110_243-492", "114-224"],  # discontinuous
    "P12931": ["88-142", "158-246", "264-524"],
    "P13489": ["7-460"],
    "P12270": [],  # no domains
    "Q0WV90": ["6-157", "347-669", "738-1081"],
}


def main():
    predictor = CCCPredictor()

    print("=" * 70)
    print("CCC Benchmark on Dev Set")
    print("=" * 70)

    for acc, expected in DEV_SET.items():
        print(f"\n{acc}")
        print(f"  Expected: {expected if expected else 'none'}")

        try:
            result = predictor.predict_from_uniprot(acc)
            predicted = [str(d) for d in result.domains]
            print(f"  Predicted: {predicted if predicted else 'none'}")
            print(f"  Method: {result.method_used}")

            # Simple comparison
            if len(expected) == len(predicted):
                print(f"  Count: MATCH ({len(expected)} domains)")
            else:
                print(f"  Count: MISMATCH (expected {len(expected)}, got {len(predicted)})")

        except Exception as e:
            print(f"  ERROR: {e}")

    print("\n" + "=" * 70)


if __name__ == "__main__":
    main()
