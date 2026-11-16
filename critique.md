# Critique of the first draft

The strongest part of my draft is that I am not only showing one cluster picture. I compare different settings, check genes that were not used to make the clusters, and keep the physical location of each spot in the analysis. The figures make the main result easy to see: a boundary can be stable in clustering while the expression change across it is gradual.

The biggest weakness is the dataset itself. It is one section from one mouse, and each spot can contain multiple cells. I cannot make a population-level claim or estimate a true cell-level transition width from this. I added the quality-matched comparison after noticing the quality associations, so it is a follow-up sensitivity check rather than something I planned before seeing the results. The step-versus-gradient comparison is simple, and some profiles may combine disconnected pieces of a boundary.

For the next version, I want to register the spots to a mouse-brain atlas or use blinded histology annotations. I also want to repeat the analysis on more 10x sections and add a spatial model that accounts for depth and mitochondrial percentage. I would analyze connected pieces of a boundary separately instead of pooling them, and I would like a more formal way to compare the step and gradient models. I also want to check the top genes against known cell-type markers before trying to interpret them.

The questions I would bring to the next class are:

- Is one tissue section okay for this first draft, or should I add multiple serial sections now?
- Do I need atlas or histology registration before calling something an anatomical boundary?
- Is the block bootstrap and held-out check enough for this class, or should I add a formal spatial regression model?
- How original does the method need to be for Option 4 if the project is mainly a careful analysis of an existing dataset?
