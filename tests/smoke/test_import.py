def test_import():
    import sandy
    assert sandy.__version__ is not None


def test_import_submodules():
    import sandy.sections
    import sandy.mcnp
    import sandy.aleph2


def test_data_directories_exist():
    from pathlib import Path
    import sandy
    base = Path(sandy.__file__).parent.resolve()
    appendix = base / "appendix"
    assert appendix.exists(), "appendix folder missing from package"

    chain = appendix / "chain_yields"
    archives = appendix / "onefile_archives"
    fycorr = appendix / "fycorr"

    assert chain.exists(), "chain_yields folder missing from wheel"
    assert archives.exists(), "onefile_archives missing from wheel"
    assert fycorr.exists(), "fycorr missing from wheel"

    # Optionally check content
    assert any(chain.glob("appendix*.txt")), "chain_yields folder does not contain 'appendix*.txt' files"
    assert any(archives.glob("*.tar.xz")), "onefile_archives should contain .tar.xz files"
    assert any(fycorr.glob("*.tar.xz")), "fycorr should contain .tar.xz files"