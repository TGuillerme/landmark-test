# Associate Editor

Comments to the Author:

This version of the manuscript improves upon the previous one, but, as evident from the two reviews plus my own, several serious problems remain.
One is that you propose a new method for testing the hypothesis that evolution is constrained by the line of least resistance, but neither reviewer, nor I, is convinced that it improves upon available methods.
That hypothesis, however, might be excessively restrictive in the case of high-dimensional data, see Hansen and Voje’s critique of the study by Berner et al. (Berner et al. 2008; Hansen and Voje 2011).
Especially when little of the variance is along PC1, the hypothesis should be extended from the “line” of least resistance to “subspaces of lesser resistance” (e.g., McGuigan et al. 2005; Hunt 2007; Martinez-Abadias et al. 2012).
Even if you choose to restrict the hypothesis of constraint to PC1, you need to demonstrate that your method works properly, such as by simulations of evolution along PC1 and at varying angles to PC1, showing that your method can distinguish the magnitude of the deviation from PC1.

> Problem 1: they are not convinced by our method. She wants us to demonstrate that the method works properly with simulation though I'm not entirely sure about how to do it from the top of my head (not following her "simulate evolution along PC1 and at varying angles to PC1"?).

The second major analytic problem is your application of it to the most divergent individuals. There is no clear rationale for inferring plasticity from the difference between PC1 of the phenotypic covariance matrix (P) and the contrast between the most divergent individuals. You present neither theory nor empirical evidence to support the claim that their discrepancy is evidence of plasticity. Like Reviewer #1, I would expect the most divergent individuals to differ in more than PC1, especially when PC1 explains so little of the variation. Another technical issue is that you seem to be treating PC1 of the whole-wombat sample (or the difference between species’ means) as if they approximate the direction of evolutionary change.    

> Problem 2: there is no link between the extreme variance (between the two extreme spec) and plasticity.
This one could be solved by either doing a rarefaction (e.g. the extremes, then then 95% then the 90%, 85%, etc...). Otherwise (or additionally) we can also do it on the other PC axis.
This bit will require some help from Emma to be sure that the generalisation to all axis is mathematically kosher.

In addition to those technical issues, all of which need to be addressed, this manuscript requires major revisions for the sake of clarity. There are many confusing passages (see the several queries by Reviewer #2, including that reviewers confusion about what the figures represent). One frequent point of confusion is your use of the term “macroevolutionary” when it is not clear whether you are talking about the actual macroevolutionary pattern (not analyzed in this study) or the expected one.  For example, you refer to “macroevolutionary inferences” being obscured by plasticity (e.g., Abstract; lines 57, 63).  Do you mean that the inference about the direction of change is obscured by the impact of plasticity on the species’ mean shape? That is how I would interpret the statement. Or do you mean that the impact of plasticity on variation within species leads to inaccurate predictions of the direction of macroevolution? That would seem to make more sense if the prediction assumes that the major axis is the major axis of G.  So it is not macroevolutionary pattern that is obscured by plasticity but rather it is the major axis of G that is obscured by plasticity.

> We might want to tone down or problems to simply talk about the patterns and not call it macroevolutionary or what not.

Similarly, at some points it is not clear whether you are talking about static allometry (the intrinsic constraint) or evolutionary allometry (the outcome of the intrinsic constraint). 

> Not sure about that one, but I guess it's an easy fix.

Your emphasis on allometry is also confusing because it is not clear why you would expect that particular constraint. You cite several studies that seem to support that expectation (because they find either that PC1 is a size axis or else they find significant allometry). However, several of them analyze size measurements and those measurements are typically correlated with each other (and with overall size) even in the absence of allometry. In studies of size measurements, size is commonly PC1 but that does not imply that allometry is a constraint.  One study does not analyze size measurements but it omits the scaling step of the Procrustes superimposition (Goswami 2007). The two analyses of shape examine geographic variation (among populations) not variation within a population. None of these studies suggests that allometry is an important constraint on marsupial evolution, but they seem to be the basis for your conclusion (line 200): “Taken together, these studies reinforce our impression that varying levels of biomechanical impacts on individual skulls can obscure the detection of macroevolutionary patterns in marsupials (and probably other mammals) with high-impact mastication.” That sounds like you are saying that allometry is not detected because of high-impact mastication (although, again, it is not clear whether it is the macroevolutionary pattern or the structure of variation, which you think is obscured by plasticity. Yet another confusing point is that you seem to be contrasting “ordination-based” hypotheses (second para, p. 24) with process-based hypotheses, but the particular hypothesis that you are testing is specifically about PC1 (of G). This is not a case of algebraic construct being treated as if it corresponds to a biological factor. Rather, it is precisely because PC1 is the main axis of variation that it is hypothesized to be a constraint. Finally, I find it confusing that you seem to treat variation in size as if it were not environmental in origin, in contrast to the effects of mechanical strain. Yet, size is obviously influenced by the environment and variation in mechanical strain can be genetic in origin. The main thesis of the classic hypothesis of mammalian craniofacial integration is that bone is shaped by effects of soft tissues, hence G is structured by those “epigenetic” interactions (Cheverud 1982). Atchley and Hall (1991) also stress that the integration of the mandible can be due to genetic covariation among muscles, which they term “epigenetic pleiotropy.”

> Reading that bit made my head hurt. I am not sure what to take out of it, maybe that we make to many guesses on the biological things going on? Though I feel this is unfair since it's kind of what they asked in the previous reviews.

Atchley, W. R. and B. K. Hall. 1991. A model for development and evolution of complex morphological structures. Biological Reviews of the Cambridge Philosophical Society 66:101-157.
Berner, D., D. C. Adams, A. C. Grandchamp, and A. P. Hendry. 2008. Natural selection drives patterns of lake-stream divergence in stickleback foraging morphology. Journal of Evolutionary Biology 21:1653-1665.
Cheverud, J. M. 1982. Phenotypic, genetic, and environmental morphological integration in the cranium. Evolution 36:499-516.
Hansen, T. F. and K. L. Voje. 2011. Deviation from the line of least resistance does not exclude genetic constraints: A comment on Berner et al. (2010). Evolution 65:1821-1822.
Hunt, G. 2007. Evolutionary divergence in directions of high phenotypic variance in the ostracode genus poseidonamicus. Evolution 61:1560-1576.
Martinez-Abadias, N., M. Esparza, T. Sjovold, R. Gonzalez-Jose, M. Santos, M. Hernandez, and C. P. Klingenberg. 2012. Pervasive genetic integration directs the evolution of human skull shape. Evolution 66:1010-1023.
McGuigan, K., S. F. Chenoweth, and M. W. Blows. 2005. Phenotypic divergence along lines of genetic variance. American Naturalist 165:32-43.


# Reviewer: 2

Comments to the Author
This is a resubmission of manuscript 18-0617, which I previously reviewed for Evolution.

The authors have put visible effort to review the manuscript according to the recommendations of the reviewers and the AE. Several of the points I had raised previously have been much improved, in particular with respect to providing a more clear focus to the manuscript, which is now more centered on a biological question / system. However, I still feel that several of the main points I had raised have not been really addressed, and I retain serious doubts about both the analytical approach taken and the biological side of the analyses. These are outlined below.

1.Regarding the micro- vs. macro-evolutionary focus of the study, and the presented objective of testing for a biomechanical constraint on morphological variation with respect to e.g. allometry,

first I do not see how such an objective can be addressed by comparing species that apparently have relatively similar diets (i.e. lines 64-65).
> This one is a fair point but is about the whole design of the study and we cannot fix that. I don't really like these kind of comments since they are not addressable in any way but by just doing another study. It's like if we'd review Beck's paper and said "ancestral states reconstructions don't work" - well they're not going to reinvent the whole field of phylogenetics...

Second, I am not sure a three-species study merits the focus of “macroevolution”.
> That's not a constructive comment.

Third, I do not see that allometry can fully be addressed – or put in the context of “constraints” – by looking at static allometry alone (i.e. across a sample of adult individuals, i.e. line 103). More philosophically, I do not see why allometric variation (or the lack of it thereof) is necessarily incompatible with a strong integration as examined here through the covariation of different structures, or with the existence of constraints.
> This is a point that they don't agree with our interpretation of the results and is maybe fair.

2.A major concern is that I am not convinced that the analytical approaches taken by the authors are robust, or are contributing to a better understanding of morphological variation.
> Both this reviewer and Zelditch use the "not convincing" argument to reject the paper. It's really poor and not constructive.

First, the argument of lines 92-96 is not self-evident to me, as other factors (e.g. sexual dimorphism, or any other unexamined factor) could be contributing to what is reflected in PC1.
> This is easy, we just add a line.

Second, it is still not clear to me why a difference between variation in PC1 vs. Procrustes shape variation is biologically interpretable the way the authors argue (same lines). PC space is a rigid rotation to the axis of major variation in Procrustes space. This makes it a mathematical, not biological, property of the data. Certainly interpretable, but I do not see how adaptation or constraint can be inferred just by comparing the two.
> I think this is the core bit that they do not like, both this reviewer and Zelditch.

Third, and most importantly, the variation ranking approach taken here still seems not to fulfill the objective presented in e.g. lines 96-99, as it seems to focus on the magnitude of variation alone, not taking variation direction into account. Methods based on trajectory analyses are available for the same purpose, so I am not sure what this new approach adds to the morphometric toolkit. Along the same lines, studying within- vs. across species covariation of shape with size and sources of morphological variation has been previously done through the common allometric component (Mitteroecker et al. 2004, JHE 46: 679-698).
> I guess we can also add these methods to the analysis to please them?

3.Again in the analytical side of things, I have a major, not previously mentioned doubt, about data analyses and what is actually presented in the manuscript. Focusing on figures 3 and 4, I get a feeling that either the figures are not correctly constructed, or not clearly explained. How can the amounts of variation be different when looking at the two extremes of the same PC axis? If I am interpreting the figure correctly, these are the theoretical shapes corresponding to the min and max of PC1, so the vectors on each landmark should have the same magnitude and exactly opposite directions. This is, of course, difficult to visualize when projecting 3D shapes on a figure, but still should be visible. Also, while I do not know much about wombat skulls, the two extremes look quite unintuitive to me with respect to the distribution of the landmarks, which makes me suspect something weird might have happened while placing the surface landmarks? This is just a guess driven by contrasting the empty area in the frontal / nasal area in one extreme and the fully landmark-populated area in the other. 
> If I understand this comment well, they basically don't get how to look at the figure. One way to do that would be to plot the landmarks for the mean specimen and then the landmarks for each max min on each side (so instead of lollypops we'll have a circle with two sticks).

4.Again on the approach proposed for inferring which areas of the skull vary the most, I still believe that the division into subsets is completely arbitrary (line 191; line 37 in p. 15), which is something the authors did not address upon my previous comment.
> I don't remember what we answered but yes. It's arbitrary I guess. But in the same way we choose to study wombats and do this analysis. We just decided this is the way we go. I don't think that arbitrarily studying the zygomatic arch is a bad thing.

5.In terms of the point I had raised previously of the way concepts and analyses are described, I still think there are numerous instances where the authors are vague or use ambiguous or confusing wording. These include, but are not limited to, lines 75-77 (why would an ecological pressure on biomechanics translate into a lack of allometry? And what type of allometry? Ontogenetic, static, evolutionary? All?); lines 87-88 (within individuals???); lines 100-101 (in what context exactly do you expect these caveats to apply? What do you mean by using within-species variation for macroevolutionary inferences?); lines 159 – 161 (what do you mean by “PC1-based patterns are significant? why would within-species PC1 be representative of a macroevolutionary trend?); line 221 (“whether the selected partitions had a significant magnitude of change”, not sure what this means).
> These seem to be another round of clarifications.

Based on these points, I am not convinced that the approach taken makes a strong inference with respect to the question of how within-species variation translates into among-species patterns, or that the analyses presented provide clear evidence on the role of allometry and biomechanical constraints in shaping morphological patterns. For this reason, I find the discussion quite overstated, with several points not really supported by data analyses (e.g. l. 150-151 p. 22; l. 154 in the same page & paragraph starting in l. 173, p. 23: concordance of PCs across species was not tested for), while others are not – as far as I can say – clearly linked to the data and analyses conducted (e.g. l. 166-168, p. 22; l. 200 – 203, p. 24).
> They are probably right for this point and the least we can do is tone down the discussion (again). However, as I was saying, it's kind of unfair when you see what other papers can get up with. I wonder what makes it so.

Reviewer: 1

Comments to the Author

This manuscript (18-0792) addresses an interesting and important question in morphological evolution: do the same factors drive micro- and macroevolutionary phenotypic variation. The parts of this manuscript comparing allometric vs biomechanical factors influencing intra- and interspecific shape variation are very good and represent an important contribution to understanding marsupial cranial evolution. However, it is hard to understand how the parts of the manuscript dealing with the permutation test for comparing PC1 variation and divergence in Procrustes space contribute to testing the evolutionary hypotheses put forth in the introduction and how they can advance the available toolkit for studying shape evolution.

If I am understanding the permutation tests correctly, you are evaluating whether the variation described by PC1 and the shape difference between the most different 2 specimens with a sample describe similar patterns of variation. You find that there is very poor correspondence between shape change on PC1 and maximal shape change in Procrustes space and conclude that within species morphospaces contain lots of “noisy” variation. The poor correspondence should not be surprising: Table 2 shows that PC1 only describes 20-25% of the variation in any of the ordinations. Within each species or genus, one would expect that the individuals with the highest pairwise Procrustes distance should be very divergent from one another in many trait dimensions (ie, they should be far apart in on many PC axes). Because PC1 describes a small amount of the overall variance, it makes sense that PC1 variance and the inter-individual differences don’t describe the same shape differences.
> What is this comment saying: just "yeah we knew that" kind of thing? I feel this does not need to be addressed.

With this in mind, it is hard for me to understand what new insights are gained by this new method. If the goal is to determine whether intraspecific and interspecific variation have similar “lines of least resistance,” perhaps something like a random skewers test would be more appropriate. If the goal is to demonstrate that PC1 scores don’t accurately summarize shape variation, this is probably an unnecessarily complex way to do it. 
> Fair enough, if that's what they want, we can drop it...

My recommendation is that either:
(1) the question being addressed by this new method needs to be more clearly stated, and why it is superior other similar methods needs to be articulated OR
(2) the sections related to the landmark variation test should be removed entirely
> We could actually follow suggestion 2: drop the landmark test and go for the classic tests suggested by the reviewer above. We can then develop the landmark variation test as a proper standalone methodological paper...

Line 9-11: “short-term biomechanical “noise” could obscure such macroevolutionary inferences”: The sentence is unclear to me. What time scale is the biomechanical signal acting on? What is the relationship between biomechanics, macroevolution, and intraspecific variation? The point is well articulated in lines 56-59, use that type of explanation in the abstract. 

Line 17-18: “does not reflect an important determinant of individual shape:” I think this is overstating the case. What is shown in the text is that there is weak correspondence between the variation between the most divergent pair of specimens and PC1. 

Line 72-85: Are the allometric and biomechanical hypotheses mutually exclusive?

Line 91: This should be “without isometric size”

Line 136: comma should be removed

Line 138: this parenthetical phrase begins with a “(“ and ends with a “,”. Choose one.

Line 228: can you provide a citation to support the idea that this is an epigenetic factor? 

Line 287: I think “wombat” should be singular

Table 5: what are the values in the cells? Area difference? Bhattacharyya coefficients?