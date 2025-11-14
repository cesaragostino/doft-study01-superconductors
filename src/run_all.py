"""Wrapper para ejecutar el master pipeline desde la carpeta src."""

from doftstudy import run_master_pipeline


def main():
    """Delegar directamente en el master pipeline."""

    run_master_pipeline.main()


if __name__ == "__main__":
    main()
