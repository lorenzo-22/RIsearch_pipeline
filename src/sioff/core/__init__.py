"""Pure core library for siOFF.

Functions here orchestrate the stateless services and return in-memory objects
(Polars DataFrames, dicts, generators). They contain no CLI concerns: no stdout
printing, no Typer/Click exceptions — they raise plain Python exceptions and write
no files. The Typer commands in ``sioff.commands`` and the public API in
``sioff.api`` are thin wrappers over this layer.
"""
