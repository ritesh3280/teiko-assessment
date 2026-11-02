"""Analysis code for the Loblaw Bio miraclib trial (Parts 2-4).

Every function takes an open ``sqlite3.Connection`` to the database built by
``load_data.py`` and returns pandas DataFrames, so the same code drives the
pipeline outputs (``run_analysis.py``) and the dashboard (``dashboard.py``).

* ``analysis.frequencies`` - Part 2: relative frequency of each population per sample.
* ``analysis.responders``  - Part 3: responders vs non-responders statistics.
* ``analysis.subsets``     - Part 4: baseline subset queries.
* ``analysis.plots``       - static figures written by the pipeline.
* ``analysis.report``      - the Markdown summary written by the pipeline.
"""
