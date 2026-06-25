"""Pure core library for RIOT.

Functions here orchestrate the stateless services and return in-memory objects
(Polars DataFrames, dicts, generators). They contain no CLI concerns: no stdout
printing, no Typer/Click exceptions — they raise plain Python exceptions and write
no files. The Typer commands in ``riot.commands`` and the public API in
``riot.api`` are thin wrappers over this layer.
"""
