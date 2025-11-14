"""Wrapper para ejecutar el master pipeline desde la carpeta src."""

try:
    from src.doftstudy import run_master_pipeline
except ModuleNotFoundError:
    # Fallback when running as ``python src/run_all.py`` without installing the package
    from pathlib import Path
    import sys

    SRC_DIR = Path(__file__).resolve().parent
    if str(SRC_DIR) not in sys.path:
        sys.path.insert(0, str(SRC_DIR))
    from doftstudy import run_master_pipeline


def main():
    """Delegar directamente en el master pipeline."""

    run_master_pipeline.main()


if __name__ == "__main__":
    main()
