#!/usr/bin/env python3
"""
Test script to verify MOFA+ installation and dependencies.

This script checks that mofapy2 is properly installed and functional.
"""

import sys


def test_import():
    """Test if mofapy2 can be imported."""
    print("Testing mofapy2 import...", end=" ")
    try:
        import mofapy2
        print("✓ OK")
        return True
    except ImportError as e:
        print("✗ FAILED")
        print(f"Error: {e}")
        print("\nTo install mofapy2, run:")
        print("  pip install mofapy2")
        print("or")
        print("  conda install -c bioconda mofapy2")
        return False


def test_entry_point():
    """Test if MOFA entry point can be imported."""
    print("Testing MOFA entry point...", end=" ")
    try:
        from mofapy2.run.entry_point import entry_point
        print("✓ OK")
        return True
    except ImportError as e:
        print("✗ FAILED")
        print(f"Error: {e}")
        return False


def test_dependencies():
    """Test required dependencies."""
    deps = ["numpy", "pandas", "h5py", "scipy"]
    all_ok = True

    print("\nTesting dependencies:")
    for dep in deps:
        print(f"  {dep}...", end=" ")
        try:
            __import__(dep)
            print("✓ OK")
        except ImportError:
            print("✗ FAILED")
            all_ok = False

    return all_ok


def test_basic_functionality():
    """Test basic MOFA functionality."""
    print("\nTesting basic MOFA functionality...", end=" ")
    try:
        import numpy as np
        from mofapy2.run.entry_point import entry_point

        # Create small dummy data
        np.random.seed(42)
        data1 = np.random.randn(10, 20)  # 10 samples, 20 features
        data2 = np.random.randn(10, 15)  # 10 samples, 15 features

        # Initialize MOFA
        ent = entry_point()
        ent.set_data_matrix(
            [data1, data2],
            views_names=["view1", "view2"],
            groups_names=["group1"],
            samples_names=[f"sample{i}" for i in range(10)],
            features_names=[
                [f"feature1_{i}" for i in range(20)],
                [f"feature2_{i}" for i in range(15)]
            ]
        )

        # Set minimal options
        ent.set_data_options(scale_groups=False, scale_views=False, center_groups=True)
        ent.set_model_options(factors=2)
        ent.set_train_options(iter=10, convergence_mode="fast", seed=42, verbose=False)

        # Build (but don't train fully)
        ent.build()

        print("✓ OK")
        return True

    except Exception as e:
        print("✗ FAILED")
        print(f"Error: {e}")
        return False


def main():
    """Run all tests."""
    print("="*60)
    print("MOFA+ Installation Test")
    print("="*60)
    print()

    tests = [
        test_import(),
        test_entry_point(),
        test_dependencies(),
        test_basic_functionality()
    ]

    print()
    print("="*60)
    if all(tests):
        print("✓ All tests passed!")
        print("MOFA+ is properly installed and ready to use.")
        print("="*60)
        return 0
    else:
        print("✗ Some tests failed.")
        print("Please install missing dependencies and try again.")
        print("="*60)
        return 1


if __name__ == "__main__":
    sys.exit(main())
