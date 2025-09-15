      seqfile = core70_XCL1_CCL20_only.phy  * pruned PHYLIP alignment (amino acids) *
     treefile = core70_XCL1_CCL20_only.phy_phyml_tree.nwk  * pruned Newick tree with branch lengths *
      outfile = codeml.out

        noisy = 3
      verbose = 1
      runmode = 0                        * user tree *

        seqtype = 2                      * 2 = amino acids *
       aaRatefile = /yezhirui/evo_probe/paml-4.10.9-linux-x86_64/dat/jones.dat            * JTT matrix *
           clock = 0
          model = 2                      * empirical AA model *
         NSsites = 0

     fix_blength = 1                     * use given branch lengths *
      cleandata = 0                      * keep gaps/missing as '?' *

     fix_alpha = 0                       * estimate gamma shape (alpha) *
          alpha = 1.0                    * initial value *
          ncatG = 4                      * JTT+Gamma(4) *

       getSE = 0
 RateAncestor = 1                        * reconstruct ancestors; output to 'rst' *
