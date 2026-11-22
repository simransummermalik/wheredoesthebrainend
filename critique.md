# Critique of the first draft

This draft has more than one cluster map: I compare different settings, check a separate set of genes, and keep each spot's location in the analysis. The result I want to look at more closely is that a boundary can stay in the same place across runs even when expression changes gradually across it.

I only used one section from one mouse, and each spot can contain several cells. That means I can't tell yet whether the same interfaces show up in other mice or how wide a transition is at the cell level. I also added the quality-matched comparison after I saw that expression change was related to depth and mitochondrial percentage, so I need to treat that result as an extra check. Some of the profiles may also combine separate pieces of the same boundary.

For the next version, I want to compare the spots with a mouse-brain atlas or histology annotations and run the analysis on more sections. I want to look at separate connected pieces of a boundary and check the top genes against known cell-type markers. I would also like to improve the comparison between the step and gradient models and account for depth and mitochondrial percentage in a spatial model.

I want to ask in class:

- Is one section enough for this draft, or should I add more now?
- Should I register the data to an atlas before I compare the results with anatomy?
- Are the held-out genes and block bootstrap enough for this project, or should I add a spatial regression model?
- For Option 4, does the research question need a new method, or is it okay to use existing methods to investigate it in a new way?
