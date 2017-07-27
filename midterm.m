% Q6: Download from Genbank a complete mitochondrial genome (mtDNA) of an animal of your
% choice. Describe the data by reporting: the complete name of the organism, accession number
% of the sequence and basic statistics including at least sequence length, base counts,
% CG-content plot. 
getgenbank('AP003425', 'sequenceOnly', true)
Error using getncbidata (line 191)
The key AP003425 was not found in the nucleotide database at this time.
Please check that the input is a valid accession number or try again.

NOTE:  This function is dependent on NCBI's Entrez tools and sequence databases. Changes to either may cause this function to
break.

Error in getgenbank (line 70)
    [varargout{1:nargout}] = getncbidata(accessnum,'fileformat','GenBank','database','nucleotide',varargin{:});
 
help bioinfo
  Bioinformatics Toolbox
  Version 4.6 (R2016a) 10-Feb-2016
 
  File I/O
    BioIndexedFile   - Read access to large text files using an index file.
    affyread         - Read Affymetrix GeneChip files. 
    affyprobeseqread - Read Affymetrix GeneChip probe sequence file.
    agferead         - Read Agilent Feature Extraction format data.
    bamindexread     - Read the index of a BAM formatted file.
    baminfo          - Return a summary of the content of a BAM file.
    bamread          - Read a BAM formatted file.
    blastread        - Read an NCBI BLAST format report file.
    blastreadlocal   - Read local BLAST format report file.
    celintensityread - Read probe intensities from Affymetrix CEL files.
    cytobandread     - Read cytogenetic banding information.
    emblread         - Read an EMBL format file.
    fastainfo        - Return a summary of the contents of a FASTA file.
    fastaread        - Read a sequence from a FASTA format file or URL.
    fastawrite       - Write a sequence to a FASTA format file.
    fastqinfo        - Return a summary of the contents of a FASTQ file.
    fastqread        - Read a FASTQ format file.
    fastqwrite       - Write sequences to a FASTQ format file.
    galread          - Read GenePix GAL file.
    genbankread      - Read a GenBank format file.
    genpeptread      - Read a GenPept format file.
    geoseriesread    - Read Gene Expression Omnibus (GEO) GSE format data.
    geosoftread      - Read Gene Expression Omnibus (GEO) SOFT format data.
    getblast         - Get a BLAST report from NCBI.
    getembl          - Get sequence data from EMBL.
    getgenbank       - Get sequence data from GenBank.
    getgenpept       - Get sequence data from GenPept.
    getgeodata       - Get Gene Expression Omnibus (GEO) data.
    gethmmalignment  - Get a multiple alignment from the PFAM database.
    gethmmprof       - Get a HMM from the PFAM database.
    gethmmtree       - Get a phylogenetic tree from the PFAM database.
    getpdb           - Get sequence data from PDB.
    gprread          - Read GenePix GPR file.
    ilmnbsread       - Read data exported from Illumina BeadStudio. 
    imageneread      - Read ImaGene format results file.
    jcampread        - Read JCAMP-DX file.
    multialignread   - Read a multiple sequence alignment file.
    multialignwrite  - Write a multiple sequence alignment file.
    mzxmlinfo        - Return information about mzXML file 
    mzxmlread        - Read mzXML file.
    pdbread          - Read a PDB format file.
    pdbwrite         - Write a PDB format file.
    pfamhmmread      - Read a PFAM format HMM profile.
    phytreeread      - Read NEWICK tree formatted file.
    saminfo          - Return a summary of the content of a SAM file.
    samread          - Read a SAM formatted file.
    soapread         - Read a SOAP aligner formatted file.
    scfread          - Read an SCF format trace file.
    sffinfo          - Return a summary of the content of an SFF file.
    sffread          - Read an SFF format file.
    sptread          - Read SPOT format file.
    tgspcinfo        - Return information about an SPC format file.
    tgspcread        - Reads a Thermo-Galactic SPC format file.
 
  Sequence Conversion
    aa2int          - Convert from amino acid to integer representation.
    aa2nt           - Convert a sequence of amino acids to nucleotides.
    dna2rna         - Convert a sequence of DNA nucleotides to RNA.
    int2aa          - Convert from integer to amino acid representation.
    int2nt          - Convert from integer to nucleotide representation.
    nt2aa           - Convert a sequence of nucleotides to amino acids.
    nt2int          - Convert from nucleotide to integer representation.
    rnaconvert      - Convert RNA structure between bracket and matrix representation.
    rna2dna         - Convert a sequence of RNA nucleotides to DNA.
    seq2regexp      - Convert a sequence that contains wildcards to a regular expression.
    seqcomplement   - Calculate the complementary strand of a DNA sequence.
    seqrcomplement  - Calculate the reverse complement of a DNA sequence.
    seqreverse      - Reverse a sequence.
 
  Sequence Statistics
    aacount         - Report amino acid counts in a sequence.
    atomiccomp      - Calculate atomic composition of a protein.
    basecount       - Report nucleotide base counts in a sequence.
    BioReadQualityStatistics - quality statistics from short-read sequences.
    codonbias       - Report codon usage per amino acid for a DNA sequence.
    codoncount      - Report codon counts in a sequence.
    cpgisland       - Locate CpG islands in a DNA sequence.
    dimercount      - Report dimer counts in a sequence.
    featurecount    - Report the number of reads mapping to genomic features.
    isoelectric     - Estimate the isoelectric point of a protein sequence.
    molweight       - Calculate molecular weight of a peptide sequence.
    nmercount       - Report n-mer counts in a sequence.
    ntdensity       - Plot nucleotide density along the sequence.
    seqwordcount    - Report word counts for a sequence.
 
  Sequence Utilities
    aminolookup      - Lookup table for peptide symbols.
    baselookup       - Lookup table for nucleotide symbols.
    cleave           - Cleave a protein with an enzyme.
    evalrasmolscript - Send Rasmol script to a molecule viewer.
    featureview      - Graphical map showing the features of a GenBank structure.
    featureparse     - Parse features from GenBank, GenPept, or, EMBL data.
    geneticcode      - Mapping for the genetic code.
    joinseq          - Join two sequences.
    ngsbrowser       - Browse short-read sequence alignment.
    molviewer        - Visualize molecules.
    oligoprop        - DNA oligonucleotide sequence properties.
    palindromes      - Find palindromes in a sequence.
    pdbdistplot      - Visualization of inter-molecular distances in a PDB file.
    proteinplot      - GUI for protein analysis.
    ramachandran     - Ramachandran plot for PDB data.
    randseq          - Generate a random sequence from a finite alphabet.
    rebasecuts       - Find restriction enzymes that cut a sequence.
    restrict         - Split a sequence at a restriction site.
    revgeneticcode   - Reverse mapping for the genetic code.
    rnaplot          - Plot RNA secondary structure.
    rnafold          - Predict secondary structure of a RNA sequence.
    seqconsensus     - Compute the consensus sequence for a set of sequences.      
    seqdisp          - Format long sequences for easy viewing.
    seqlogo          - Display sequence logos for DNA and protein sequences.
    seqmatch         - Find matches for every string in a library.
    seqprofile       - Compute the sequence profile of a multiple alignment.
    seqshoworfs      - Graphical display of Open Reading Frames in a sequence.
    seqshowwords     - Graphical display of words in a sequence.
    seqviewer        - Visualize biological sequences.
 
  Sequence Alignment
    align2cigar      - Calculate the CIGAR of aligned sequences.
    blastformat      - Run local version of BLAST formatdb.
    blastlocal       - Run local version of BLAST.
    blastncbi        - Generate a remote NCBI BLAST request.
    bowtie           - Map short reads using the Burrows-Wheeler transform.
    bowtiebuild      - Generate an index using Burrows-Wheeler transform.
    cigar2align      - Align sequences using a CIGAR string.
    localalign       - Find optimal and suboptimal local alignments of two sequences.
    multialign       - Progressive multiple sequence alignment.
    ngsbrowser       - Interactive browser for exploring short read sequence alignment data.
    nwalign          - Needleman-Wunsch global alignment.
    profalign        - Needleman-Wunsch global alignment of two profiles.
    seqalignviewer   - Visualize and edit a sequence alignment.
    seqdotplot       - Create a dotplot of two sequences.
    showalignment    - Visualization of pairwise sequence alignment.
    swalign          - Smith-Waterman local alignment.
 
  Sequence Data Containers
    BioRead          - Class representing a collection of sequences with quality scores.
    BioMap           - Class representing a collection of sequences with alignment information.
    GFFAnnotation    – Class representing a collection of GFF annotations.
    GTFAnnotation    – Class representing a collection of GTF annotations.
 
  Statistical Learning
    classify         - Discriminant analysis (Statistics and Machine Learning Toolbox).
    classperf        - Performance of a classifier with diagnostic tests.
    classregtree     - Classification tree fitting (Statistics and Machine Learning Toolbox).
    crossvalind      - Cross-validation index generation.
    fitcknn          - Fit a K-Nearest neighbor classifier (Statistics and Machine Learning Toolbox).
    fitcsvm          - Fit a Support Vector Machine classifier (Statistics and Machine Learning Toolbox).
    kmeans           - K-means clustering (Statistics and Machine Learning Toolbox).
    kmedoids         - K-medeoids clustering (Statistics and Machine Learning Toolbox).
    knnimpute        - Impute missing data using the nearest neighbor method.
    metafeatures     - Attractor metagene learning.
    nbintest         - Unpaired hypothesis test for count data with small samples.
    optimalleaforder - Reorders a hierarchical binary cluster tree. (Statistics and Machine Learning Toolbox)
    randfeatures     - Randomized subset feature selection.
    rankfeatures     - Ranks key features by class separability criteria.
    samplealign      - Aligns two data sets containing sequential observations. 
 
  Protein Analysis
    aacount          - Show the amino acid composition of a protein sequence.
    aminolookup      - Lookup table for peptide symbols.
    atomiccomp       - Calculate atomic composition of a protein.
    cleave           - Cleave a protein with an enzyme.
    cleavelookup     - Display cleavage rules of enzymes or compounds.
    evalrasmolscript - Send Rasmol script to a molecule viewer.
    isoelectric      - Estimate the isoelectric point of a protein sequence.
    molviewer        - Visualize molecules.
    molweight        - Calculate molecular weight of a peptide sequence.
    pdbdistplot      - Visualization of inter-molecular distances in a PDB file.
    pdbsuperpose     - Superpose the 3-D structures of two proteins.
    pdbtransform     - Apply a linear transformation to the 3D structure of a molecule.
    proteinpropplot  - Plot hydrophobicity and other properties of a sequence.
    proteinplot      - GUI for protein analysis.
    ramachandran     - Ramachandran plot for PDB data.
 
  Trace tools
    scfread         - Read SCF format trace data.
    traceplot       - View nucleotide trace plots.
 
  Profile Hidden Markov Models
    gethmmalignment - Get a multiple alignment from the PFAM database.
    gethmmprof      - Get a HMM from the PFAM database.
    gethmmtree      - Get a phylogenetic tree from the PFAM database.
    hmmprofalign    - Sequence alignment to a profile HMM.
    hmmprofestimate - Estimate the parameters of a profile HMM.
    hmmprofgenerate - Generate a random sequence from a profile HMM.
    hmmprofmerge    - Align the output strings of several profile alignments.
    hmmprofstruct   - Create a profile HMM structure.
    pfamhmmread     - Read a PFAM format HMM profile.
    showhmmprof     - Plot an HMM profile.
 
  Phylogenetic Tree Tools
    dnds            - Estimate synonymous and nonsynonymous substitution rates.
    dndsml          - DNDS using maximum likelihood.
    phytree         - Class representing a phylogenetic tree object.
    phytreeread     - Read NEWICK tree formatted file.
    phytreeviewer   - Visualize and edit phylogenetic trees.
    phytreewrite    - Save a phylogenetic tree object as a NEWICK format file.
    seqlinkage      - Construct a phylogenetic tree from pairwise distances.
    seqneighjoin    - Neighbor-joining for phylogenetic tree reconstruction.
    seqpdist        - Pairwise distance between sequences.
 
  Phylogenetic Tree Methods
    phytree/cluster      - Construct clusters from a phylogenetic tree.
    phytree/get          - Get information about a phylogenetic tree object.
    phytree/getbyname    - Select branches and leaves by name.
    phytree/getcanonical - Calculates the canonical form of a phylogenetic tree.
    phytree/getmatrix    - Converts a tree into a relationship matrix.
    phytree/getnewickstr - Creates a NEWICK formatted string.
    phytree/pdist        - Compute the pairwise patristic distance.
    phytree/plot         - Render a phylogenetic tree.
    phytree/prune        - Reduce a phylogenetic tree by removing branch nodes.
    phytree/reroot       - Changes the root of a phylogenetic tree.
    phytree/reorder      - Changes the leaf order of a phylogenetic tree.
    phytree/select       - Select tree leaves and branches.
    phytree/subtree      - Extracts a subtree.
    phytree/view         - View a phylogenetic tree in phytreeviewer.
    phytree/weights      - Tree-based sequence weights.
 
  Microarray Data Containers
    bioma.ExpressionSet     - Class to contain microarray gene expression experiment data.
    bioma.data.DataMatrix   - Two dimensional data array with row and column names.
    bioma.data.ExptData     - Data container class to store experiment data values.
    bioma.data.MetaData     - Metadata collection of variables and their description.
    bioma.data.MIAME        - Class for storing a description about a microarray experiment.
 
  Microarray Data Analysis and Visualization
    cghcbs          - Compute circular binary segmentation on array CGH data.
    cghfreqplot     - Frequency plot of copy number alterations.
    chromosomeplot  - Plot chromosome ideograms.
    clustergram     - Clustergram plot.
    HeatMap         - False color image of numeric array.
    maboxplot       - Box plot of microarray data.
    mafdr           - Compute false discovery rates of gene expression data.
    maimage         - Pseudocolor plot of microarray spatial data.
    mairplot        - Intensity plot of microarray signals.
    maloglog        - Log-log plot of microarray data.
    mapcaplot       - Principal Component plot of expression profile data.
    mattest         - Unpaired student's t-test of microarray data.
    mavolcanoplot   - Volcano plot of expression profile data.
    microplateplot  - Creates a visualization of a microtiter plate.
    redbluecmap     - Generate red and blue colormap.
    redgreencmap    - Generate red and green colormap.
 
  Microarray Normalization and Filtering
    affygcrma         - Perform GCRMA procedure for multiple Affymetrix GeneChips.
    affyinvarsetnorm  - Invariant set normalization of Affymetrix probe-level data.
    affyprobeaffinities - Compute probe affinities of Affymetrix GeneChip.
    affyrma           - Perform RMA procedure for multiple Affymetrix GeneChips.
    exprprofrange     - Calculate range of expression profiles.
    exprprofvar       - Calculate variance of expression profiles.
    gcrma             - GCRMA measure gene expression of Affymetrix microarray data.
    gcrmabackadj      - GCRMA background adjustment of Affymetrix probe-level data.
    geneentropyfilter - Remove genes with entropy expression values.
    genelowvalfilter  - Remove genes with low expression values.
    generangefilter   - Remove genes with small expression ranges.
    genevarfilter     - Remove genes with small expression variance.
    mainvarsetnorm    - Rank invariant set normalization.
    malowess          - Lowess normalization.
    manorm            - Normalization by scaling and centering.
    quantilenorm      - Quantile normalization.
    rmabackadj        - RMA background adjustment of Affymetrix probe-level data.
    rmasummary        - RMA summarization of multiple Affymetrix microarray data.
    zonebackadj       - Zone based background adjustment of Affymetrix probe-level data.
 
  Microarray Utility Functions
    affysnpintensitysplit - Splits Affymetrix SNP probe intensity matrix for allele A and B.
    affysnpquartets       - Create a table of SNP quartets for a probe set.
    ilmnbslookup          - Look up Illumina probe sequence and annotation information.
    magetfield            - Extract data from microarray structure.
    probelibraryinfo      - Get library information for a probe.
    probesetlink          - Show probe set information from NetAffx.
    probesetlookup        - Get gene information for a probe set.
    probesetplot          - Plot probe set values.
    probesetvalues        - Get probe set values from CEL and CDF information.
 
  Gene Ontology Functions and Methods
    geneont                        - Creates a Gene Ontology (GO) object
    geneont.geneont/getancestors   - Finds the ancestors of a GO term
    geneont.geneont/getdescendants - Finds the descendents of a GO term
    geneont.geneont/getmatrix      - Converts a GO Object into a relationship matrix
    geneont.geneont/getrelatives   - Finds the related terms for a GO term
    goannotread                    - Extract data from microarray structure.
    num2goid                       - Converts numeric values to GO IDs
 
  Bioanalytics and Mass-Spectrometry Preprocessing and Visualization
    isotopicdist      - Calculate isotope mass distribution and density function.
    jcampread         - Read JCAMP-DX file.
    msalign           - Signal calibration and alignment by reference peaks.
    msbackadj         - Background estimation and correction.
    msdotplot         - Create a dot plot of an LCMS or GCMS dataset.
    msheatmap         - Heat map image of a set of spectra.
    mslowess          - Non-parametric smoothing using Lowess method.
    msnorm            - Normalization of a set of spectra.
    mspalign          - Peak binning and dynamic programming peak alignment.
    mspeaks           - Peak detection with wavelet denoising.
    msppresample      - Signal resampling from peak information.
    msresample        - Resample with antialias filtering.
    mssgolay          - Least-squares polynomial smoothing.
    msviewer          - Plot a spectrum or a set of spectra.
    mzcdf2peaks       - Convert an mzCDF structure to a list of peaks.
    mzcdfinfo         - Return information about netCDF file.
    mzcdfread         - Read netCDF file with mass-spectrometric data.
    mzxml2peaks       - Convert an mzXML structure to a list of peaks.
    mzxmlinfo         - Return information about mzXML file.
    mzxmlread         - Read mzXML file.
    samplealign       - Constrained dynamic programming alignment and warping.
 
  Graph Theory Algorithms
    graphallshortestpaths - Find distance of all shortest paths.
    graphconncomp         - Strong and weak connected components.
    graphisdag            - Check if graph is DAG.
    graphisomorphism      - Map between two isomorphic graphs.
    graphisspantree       - Check if graph is a spanning tree.
    graphmaxflow          - Max-flow (and min-cut) algorithm.
    graphminspantree      - Find the minimal spanning tree.
    graphpred2path        - Covert from a predecessor list to a path.
    graphshortestpath     - Find the shortest path between two nodes.
    graphtopoorder        - Topological order of a DAG.
    graphtraverse         - Depth first search and breadth first search.
  
  Graph Visualization Methods
    biograph                           - Create a bioinformatics graph object.
    biograph.biograph/dolayout         - Calculate node and edge positions.
    biograph.biograph/getmatrix        - Get the relationship matrix.
    biograph.biograph/getnodesbyid     - Get handles to nodes.
    biograph.biograph/getedgesbynodeid - Get handles to edges.
    biograph.biograph/view             - Render a graph in its viewer.
    biograph.node/getancestors         - Find ancestors.
    biograph.node/getdescendants       - Find descendants.
    biograph.node/getrelatives         - Find neighbors.
   
  Scoring Matrices
    blosum            - BLOSUM family of matrices.
    dayhoff           - Dayhoff matrix.
    gonnet            - Gonnet variation on PAM250.
    nuc44             - Nuc44 nucleotide matrix.
    pam               - PAM family of matrices.
 
  Tutorials, demos and examples.
    acghhmmdemo        - Bayesian hidden Markov modeling of array CGH data.
    affydemo           - Example of working with Affymetrix GeneChip data.
    affypreprocessdemo - Preprocessing Affymetrix microarray data at the probe level.
    affysnpcnvdemo     - Analyzing Affymetrix SNP arrays for DNA copy number variants.
    aligndemo          - Basic sequence alignment tutorial demo. 
    alignscoringdemo   - Tutorial showing the use of scoring matrices. 
    alignsigdemo       - Demo of how to estimate the significance of alignments.
    bacacghdemo        - Detecting DNA copy number alteration in array-based CGH data.
    birdfludemo        - Investigating the Bird Flu virus (H5N1).
    biodbdemo          - Example of connecting to local databases.
    biodistcompdemo    - Batch processing through sequential and parallel computing.
    biographdemo       - Working with BIOGRAPH objects.
    biomemorymapdemo   - Using memory mapping to work with whole genome data.
    bioperldemo        - Example of calling Bioperl functions.
    cancerdetectdemo   - Data mining analysis for mass spectrometry profiles.
    chipseqpedemo      - Exploring protein-DNA binding sites from paired-end ChIP-seq data.
    clustergramdemo    - Clustergram functionality examples.
    cnsgeneexpdemo     - Exploring gene expression data.
    diffprotdemo       - Differential proteomics and metabolomics analysis of LCMS.
    evemotifdemo       - Identifying over-represented regulatory motifs. 
    dndsdemo           - Analyzing synonymous and non-synonymous substitution rates.
    geneontologydemo   - Example of working with Gene Ontology data.
    graphtheorydemo    - Example of working with graph theoretic functions.
    gsedemo            - Example of working with Gene Expression Omnibus Series data.
    gutmicrobiomedemo  - Analyzing the human distal gut microbiome.
    hmmprofdemo        - HMM profile alignment tutorial example.
    hivdemo            - Analyzing the origin of the HIV with phylogenetic trees.
    illuminagedemo     - Example of working with Illumina BeadChip data.
    ilmnsolexademo     - Working with Illumina/Solexa next-generation sequencing data.
    lcmsdemo           - Visualizing and preprocessing LCMS data.
    maexptdemo         - Examples of working with microarray experiment data structures.
    mbdseqdemo         - Exploring genome-wide differences in DNA methylation profiles.
    metagenomicdemo    - Metagenomic analysis of a Sargasso Sea Sample.
    molviewerdemo      - Visualizing the three-dimensional structure of a molecule.
    mousedemo          - Microarray normalization and visualization example.
    msgademo           - Mass spectra data analysis with Genetic Algorithms.
    mspreprodemo       - Preprocessing of raw mass spectrometry data.
    ncbieutilsdemo     - Accessing NCBI Entrez Databases with E-Utilities.
    phybootdistdemo    - Confidence estimation of trees by bootstrapping.
    primerdemo         - Primer design tutorial example.
    primatesdemo       - Building a phylogenetic tree for the hominidae species.
    rnademo            - Predicting and visualizing the secondary structure of RNA.
    rnaseqdedemo       - Identifying differentially expressed genes from RNA-seq data.
    sarsdemo           - Reconstructing the origin and the diffusion of the SARS epidemic.
    secstructnnetdemo  - Predicting protein secondary structure using a neural network.
    seqstatsdemo       - Sequence statistics tutorial example.
    wholegenomedemo    - Comparing whole genomes.
    yeastdemo          - Microarray data analysis example.
 

fastaread('~/Downloads/hippoAmphibius.fasta')

ans = 

      Header: 'AP003425.1 Hippopotamus amphibius mitochondrial DNA, complete genome'
    Sequence: 'GTTAACGTAGCTCAAACACCCAAAGCGAGGCACTGAAAATGCCTAGATGGGCTCACCCAGCCCCGTAAACAGACAGGTTTGGTCCCAGCCTTTCTGTTAATTTTTA?'

hippo = ans.Sequence

hippo =

GTTAACGTAGCTCAAACACCCAAAGCGAGGCACTGAAAATGCCTAGATGGGCTCACCCAGCCCCGTAAACAGACAGGTTTGGTCCCAGCCTTTCTGTTAATTTTTAATAGAATTACACATGCAAGTATCCACACCCCAGTGAGAATGCCCTCTAAATCACCCGGATCAAAAGGAGCGGGTATCAAGCACACCATACACCGTAGCTCAAAACGCCCTGCTCAGCCACACCCCCACGGGAAACAGCAGTGACCAAAATTAAGCCATGAACGAAAGTTTGACTAAGCCATATTAACCAGAGTTGGTAAATCTCGTGCCAGCCACCGCGGTCATACGATTGACTCAAACTAACAGAAATACGGCGTAAAGCGTGTTAAAGAATTAAAAGTACAAATAAAGTTAAATTCTAACTAAACTGTAAAAAGCCCTAACTAGAATAAAAATCTACTACAAAAGTGACTTTAACAATACTGACCACACGATAGCTAAGACCCAAACTGGGATTAGATACCCCACTATGCTTAGCCCTAAACACAGATAATTCCAAAAACAAAACTATTCGCCAGAGTACTACTAGCAACAGCTTAAAACTCAAAGGACTTGGCGGTGCTTCATACCCCTCTAGAGGAGCCTGTTCTATAATCGATAAACCCCGATAAACCTCACCAACCCTTGCTAATCCAGTCTATATACCGCCATCTCCAGCAAACCCTAAAAAGGACTAAAAGTAAGCTCAACTATTACACATAAAGACGTTAGGTCAAGGTGTAACCTATGGGCTGGGAAGAAATGGGCTACATTTTCTAGAACAAGAACACAACCCACCCGAACGAAAACTCCTATGAAAGCTAGGAACTAAAGGAGGATTTAGTAGTAAATCAAGAGTAGAGTGCTTGATTGAACAAGGCCATGAAGCACGCACACACCGCCCGTCACCCTCCTCAAATAACACAAACCCATAAACATAATAGCTAGACAAAACCAAATGAGAGGAGACAAGTCGTAACAAGGTAAGCATACTGGAAAGTGTGCTTGGACAAATCAAAGCATAGCTTAGATTTAAAGCATCTAGTTTACACCCAGAAGATTTCACAATAAGTGAATGCCTTGAACTAAAGCTAGCCCAACTACCCCACCCCACTATAAATAAAACAAAGCATTTAATCACCAACTTAAAGTATAGGAGATAGAAATCCGAGACCATTTGGCGCAATAGAGATAGTACCGTAAGGGAATGATGAAAGAAGCAACAAAAGTATCAAAAAGCAAAGATCACACCTTGTACCTTTCGCATAATGATTTAACTAGCAAAACCTTAGCAAAGAGAACTTAAGCTAAGCCACCCGAAACCAGACGAGCTACTTATGAACAGTTTACCAGAACAAACTCATCTATGTGGCAAAATAGTGAAAAGATTTATAAGTAGAGGTGAAAAGCCTAACGAGCCTGGTGATAGCTGGTTGTCCAGAAAAAGAATCTAAGTTCAATTTCAAAATTACCAAAAGCCCAAATCAAGCCTAATGTAATTTTAAAATTTAGTCTAAAAAGGTACAGCTTTTTAGACAAAGGATACAACCTTGACTAGTGAGTAAATTCAAACAACACCATAGTTGGCCTAAAAGCAGCCACCAATTAAGATAGCGTTCAAGCTCGACAATAAGCACAACCTAATTCCCATATCAAGCAAACAACTCCTAGCCTAACTACTGGACTATTCTATTCAAACATAGAAGCAATACTGTTAATATGAGTAACAAGAAACAATTCTCCCAGCACAAGCCTACGTCAGCAACTGATAATATACTGACAGTTAACAACAAGTAAATATAACCTAACACTAAAACATTTACCCCATACACTGTTAACCCAACACAGGCGTGCACCCAAGGAAAGATTAAAAAAAGCAAAAGGAACTCGGCAAACACAAACCCCGCCTGTTTACCAAAAACATCACCTCTAGCATCACTAGTATTAGAGGCACTGCCTGCCCAGTGACAACTGTTAAACGGCCGCGGTATCCTGACCGTGCAAAGGTAGCATAATCATTTGTTCCTTAAATAGGGACTTGTATGAATGGCCACACGAGGGTTTTACTGTCTCTTGCTTTCAATCAGTGAAATTGACCTTCCCGTGAAGAGGCGGGAATAATACAATAAGACGAGAAGACCCTGTGGAGCTTCAATTAGCTGACTCAATAAAAACAAAATAAACCCGCAAGGCACAATAAAATCCTATATGAGTCAGTAATTTTGGTTGGGGTGACCTCGGAGAAAAAAGAATCCTCCGAGTGATAAAAATCTAGACTCACCAGTCAAAACATAACAACACTCATTGACCCAAAACCTTTGATCAACGGAACAAGTTACCCCAGGGATAACAGCGCAATCCTATTCTAGAGTCCATATCGACAATAGGGTTTACGACCTCGATGTTGGATCAGGACATCCCAATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAAAGTCCTACGTGATCTGAGTTCAGACCGGAGCAATCCAGGTCAGTTTCTATCTATTATACATTTCTCCCAGTACGAAAGGACAAGAGAAATGAGGCCTACTTCCAATAAGCGCCTTAAAACTAATTAATGATATAGTCTTAACTTAATTAAATAGTATAAATATACCAGCCCTAGACCAGGGCACAGTTGCGATGGCAGAGCCCGGTAATTGCATAAAACTTAAGCCTTTACACCAGAGGTTCAAATCCTCTTCACAACAAAATGTTTATCATCAACACCCTTATACTTGTCGCACCCATCCTTCTAGCCATAGCATTCCTAACACTAGTCGAACGAAAAATCCTAGGGTATATACAACTCCGAAAAGGCCCAAACATTGTAGGACCATACGGCCTACTACAACCCTTCGCCGATGCAATCAAGCTCTTCACAAAAGAACCCCTACGACCATCCACCTCCTCTGTCTCCATATTCATCATCGCACCAATCCTAGCTCTAACTCTAGCCCTAACAATATGAATCCCCCTACCCATACCATACCCCCTCATCAACATAAACTTGGGCGTACTATTTATACTAGCCATATCTAGCCTAGCTGTATATTCCATCCTATGATCCGGATGAGCCTCCAACTCGAAATATGCATTAATTGGCGCCCTACGAGCAGTAGCACAAACAATCTCATATGAAGTAACCCTAGCAATCATCCTCCTGTCCATTCTCCTAATAAACGGGTCCTTCACATTATCAACCCTCATCACAACACAGGAAAAACTGTGACTAATCTTCCCTTCATGACCACTGGCCATAATATGATTCATCTCGACCCTAGCAGAGACTAACCGAGCCCCATTCGACCTCACAGAAGGAGAATCCGAACTTGTATCAGGCTTCAACGTAGAGTATGCAGCAGGGCCATTCGCTATATTCTTCATAGCAGAATACATCAACATCATCATAATAAATGCCTTCACAACAGTCCTATTCCTAGGTGCATACCACAACCCATACTTACCAGAACTCTACACAATCAACTTTACCATCAAAACACTACTACTAACAATATCCTTCCTGTGAATCCGAGCATCCTACCCACGATTCCGATACGACCAACTAATGCACCTTCTATGAAAAAGCTTCCTACCATTAACACTAGCCCTGTGCATATGACACGTATCACTCCCCATTATAACATCAAGCATCCCCCCTCAAACATAAGAAATATGTCTGACAAAAGAATTACTTTGATAGAGTAAATAATAGAGGTTCAAGCCCTCTTATTTCTAGAACTATAGGAATCGAACCTACTCCTGAGAACTCAAAATTCCCTGTGCTACCAATCACACCCCGTTCCACAGTAAGGTCAGCTAAACAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTTAAACCCCTCCCGTACTAATAAATCCCTTCGTCTCTGTCATTATCTACACAACCATCATCCTGGGGACCATGATTGTAATAACTAGCTCTCACTGACTATTAACTTGAACCGGGTTCGAAATAAACATGCTAGCTATTATCCCCATCATAATGAAATCACCCAACCCACGAGCCACAGAAGCCTCCGTTAAATACTTCATAACCCAAGCCACCGCCTCCATGCTACTCATACTAGCAGTCATTATCAATCTACTATACTCCGGACAATGAACAGTCATAAAAATACTTAATCCAACAGCATCTATAATTATAACGATAGCTCTTGCCATAAAATTAGGACTATCCCCCTTCCACTTCTGAGTCCCTGAAGTAACACAAGGCATCCCCTTAACAGCAGGCTTAATCTTACTGACATGACAAAAACTCGCACCCCTATCCATTCTCTACCAAATCTCCCCATCAATCAATCCAAACCTAATCCTAACTATATCAATGCTATCCATCCTAGTCGGCGGATGAGGCGGGCTAAACCAAACTCAATTACGAAAAATCATAGCGTACTCATCTATTGCCCATATGGGGTGGATAGCAGCCATCCTAATCTACAACCCAACCATAACCATCCTAAACCTAACAATCTACCTCATAACAACCTTCACAATGTTTACAATATTTGCACTCAACTCAACCACCACTACCCTTTCCCTATCACACACATGAAACAAAACCCCCATTATTACAACCCTTATACTCACTATTCTACTATCAATAGGAGGACTACCCCCACTAACAGGCTTCGTACCAAAATGAATAATCATCCAAGAAATAACAAAAAACGATAGCATTATCCTACCCACACTAATAGCCATCATAGCACTCCCCAACCTATATTTCTACATACGACTCACCTACTCTACAGCACTAACCATATTTCCTTCCTCAAACAACATAAAAATAAAATGACAATTTGAAGCCTCAAAACACAAAACACTCCTGCCAACAATAATCATCCTCTCCACCATACTCCTTCCCCTCACACCAATACTAGTAGTACTAGACTAGAGGTTTAGGTTACCTAGACCAAGAGCCTTCAAAGCTCCAAGCAAGTATATATTACTTAATCTCTGCTCAATAGAGACTGCAAGACCGCACCTCACATCAATTGAATGCAACTCAACTGCTTTTATTAAGCTAAGTCCCTACTAGATTGGTGGGGCATGCTCCCCACGAACTTTTAGTTAACAGCTAAATACCCTAATCAACTGGCTTCAATCTACTTCTCCCGCCGCGGGGAAATAAAGGCGGGAGAAGCCCCGGCAGAATTGAAGCTGCTTCTTTGAATTTGCAATTCAATATGAATATTTCACTACAGGACCTGGCAAAAAGAGGACTCAACCCCTGTACTTAGATTTACAGTCTAATGCTTACTCAGCCATTTTACCCATGTTCATAAACCGCTGACTATTCTCAACCAACCACAAAGACATCGGTACACTATATCTACTATTCGGCGCCTGAGCTGGCATAGCAGGCACTGGCCTGAGCCTACTAATCCGCGCCGAACTGGGTCAACCTGGCACACTATTAGGAGATGACCAAATTTACAACGTAGTTGTCACAGCCCACGCATTTGTAATAATTTTCTTTATAGTTATACCAATTATGATTGGCGGGTTCGGAAACTGACTTGTTCCACTAATAATCGGAGCCCCTGATATGGCCTTTCCTCGAATAAATAACATAAGCTTCTGACTACTCCCTCCCTCCTTCCTACTACTATTAGCATCCTCCATGGTAGAAGCAGGGGCAGGAACAGGTTGGACCGTCTATCCCCCTTTAGCCGGAAATTTAGCCCATGCTGGAGCCTCTGTAGATCTTACAATTTTCTCCCTCCACCTAGCCGGAGTCTCCTCTATTCTAGGTGCAATCAACTTCATTACCACCATCATCAACATGAAACCACCCGCTATATCTCAGTATCAAACCCCACTGTTTGTCTGATCAGTCCTAATCACGGCTGTGCTACTTCTACTCTCCCTACCTGTTTTAGCAGCAGGTATTACTATGCTACTCACAGATCGAAACCTAAATACCACCTTCTTTGACCCTGCAGGAGGAGGCGACCCTGTCCTTTATCAACACCTATTCTGATTCTTCGGACACCCCGAAGTATACATCCTGATCCTCCCCGGCTTCGGAATAATCTCGCACATTGTAACATACTACTCCGGAAAAAAAGAACCTTTTGGGTACATAGGCATAGTCTGAGCTATAATATCCATCGGGTTCCTAGGATTTATTGTATGAGCCCATCACATATTTACAGTAGGTATAGACGTCGACACCCGAGCATACTTCACATCCGCCACTATAATTATCGCCATCCCCACAGGAGTAAAAGTATTCAGCTGACTAGCAACACTGCATGGAGGGAACATCAAATGGTCCCCTGCTATGATGTGAGCCCTAGGCTTTATTTTCCTATTCACAGTGGGTGGCCTAACAGGTATTGTTTTAGCCAACTCATCCCTAGACATCGTCCTTCACGACACCTATTACGTAGTAGCCCATTTCCATTACGTGCTCTCAATAGGCGCCGTCTTCGCTATCATAGGAGGCTTCGTACACTGATTCCCACTATTTTCAGGATATACACTCAATGACACATGAGCAAAAATCCACTTCGTAATCATGTTCGTGGGGGTCAATCTAACTTTCTTCCCACAGCATTTCTTAGGCCTATCCGGAATGCCCCGACGATACTCCGACTACCCAGACGCCTATACAACATGAAACACTATCTCCTCAATAGGTTCTTTTATCTCACTAACAGCTGTAGTACTAATAGTGTTCATCATTTGAGAGGCATTTGTCTCCAAACGAGAAGTCTTGGCTGTAGATCTAACTACAACCAACTTAGAGTGACTAAACGGGTGCCCTCCACCATACCACACATTTGAAGAACCCGCATACGTGAACTTAACTAGCCAAAACAAGAGAGGAAGGAATCGAACCTCCTCCTGTTGGTTTCAAGCCAACATCATAACCACTATGTCTCTCTCCATAAACGAGGTATTAGTAAAAATTACATAACTTCGTCAAAGTTAAGTTACAGGTGAAAACCCTGTCTACCTCCATGGCATATCCCCTCCAACTAGGCTTTCAAGATGCAGTATCACCCATTATAGAAGAACTACTGTATTTTCACGACCACACGCTAATAATCGTATTCCTAATCAGCTCACTAGTCCTTTACATTATTACACTAATACTGACTACCAAACTAACCCACACAAACACCATAAATGCACAAGAGGTAGAAACTGTCTGAACAATCCTACCAGCCATTATCCTTATCCTAATTGCACTGCCATCTCTGCGAATCCTCTATATAATAGACGAAATTAACAACCCCTCCCTGACCGTAAAAACTATGGGCCACCAATGATACTGAAGTTACGAGTATACAGATTATGAAGACCTAAACTTTGACTCCTACATAGTCCCAACATCAGACCTAAAGCCGGGGGACCTACGACTCCTAGAAGTAGATAACCGAGTCGTCCTACCCATAGATGTAACAGTTCGAATACTAATCTCATCAGAAGACGTACTACACTCCTGAGCCGTGCCATCACTAGGTCTAAAAACAGATGCCATTCCAGGACGATTGAACCAAACAACCTTAATATCAACACGACCCGGACTATTTTACGGACAGTGCTCCGAAATCTGTGGTTCCAACCACAGCTTTATGCCCATTGTCCTAGAATTAGTCCCACTGCAAACTTTCGAAAAATGAACCGCATCCCTATTATAGACTCATTAAGAAGCTAGTAGCGCTAACCTTTTAAGTTAGAGACTGAGAGCCAAGCTCTCCTTAATGACATGCCACAACTAGACACATCAACATGATTTACCACCATCCTATCCATATTTCTGACCCTATTTATTATCTTTCAACTGAAAATCTCAAAACACACCTGCCACCCAAACCCTGAGACTACTCTTCCCATAACACAAAAACAGCCTACCCCCTGAGAAACGAAATGAACGAAAATCTATTCGCCTCTTTCATTACCCCTACAATCCTAGGCCTACCCCTAGTCACCCTAATCATTATGTTCCCAAGCATACTATTCCCAGCACCCACCCGTCTAATTACTAATCGCTTAGTCTCCATTCAACAATGACTAATCCAGCTCGTATCAAAACAAATAATGAACATCCACAACCACAAAGGACAAACTTGAACACTAATATTAATATCCCTTATCCTATTCATCGGCTCAACAAACCTCCTGGGACTCCTGCCACACTCATTCACACCCACCACACAACTCTCAATAAACTTAGGCATAGCCATCCCCCTGTGAGCAGGCACTGTAATCATAGGCTTCCGTAACAAAACAAAAATCTCCCTGGCCCACTTTTTACCCCAAGGAACACCCACACCCCTAATCCCCATGCTAGTAATCATTGAGACAATCAGCCTATTTATCCAACCGATAGCACTAGCCGTACGACTAACGGCAAACATCACGGCAGGACACCTACTAATGCACCTAATCGGAGGAGCAACCCTCGCATTAATAAACATCAGCATAACCACCGCCCTTATCACGTTCATCATCCTAGTCTTACTAACAGCTCTAGAGTTTGCCGTTGCCATAATCCAAGCGTACGTCTTCACCCTACTAGTAAGTCTATACTTACACGATAACACATAATGACCCACCAAACCCACGCATACCATATAGTAAACCCAAGTCCCTGACCTCTTACAGGAGCCCTCTCAGCCCTACTAATAACGTCGGGCCTAACCATATGATTCCACTTTAACTCCCTTATCCTACTGACGACAGGACTAGTTACCAATATCCTAACAATATATCAGTGATGACGAGATGTAATCCGAGAAAGCACCTTTCAAGGCCACCACACACCAGTCGTACAAAAAGGACTTCGCTACGGAATAGTCCTATTTATTATCTCCGAAGTCCTATTTTTCACAGGCTTCTTCTGAGCCTTTTACCACTCAAGCCTCGCTCCTACTCCTGAACTAGGCGGATGTTGACCACCCACAGGCATCAACCCTTTAAACCCACTAGAGGTGCCACTTCTAAACACCTCCGTTCTACTAGCCTCTGGTGTCTCCATTACCTGAGCCCACCACAGTTTAATAGAAGGCAATCGAAAACAAATACTCCAAGCCCTCTTCATCACAATCGCCCTAGGTGTGTACTTCACACTACTGCAAGCTTCAGAATACCATGAAGCCTCCTTTACAATCTCAGATGGGGTTTACGGCTCAACTTTCTTTGTAGCCACAGGCTTTCATGGATTACATGTAATTATTGGCTCCACTTTCCTAATTGTATGCTTCCTACGCCAACTAAAATTCCACTTCACGTCAGATCACCACTTTGGCTTCGAGGCCGCCGCCTGATACTGACACTTCGTAGATGTAGTCTGACTATTCCTCTACGTGTCCATTTATTGATGAGGTTCATAGTTCTTTTAGTATTAAAACAGTACAGCTGACTTCCAATCAGCTAGCCTCAGTGCACTCTGGGAAGGAACAATTAACCTAATAATAGCACTACTAACAAACACCGCACTAGCCTCTTTATTGGTCCTCATTGCCTTCTGACTCCCACAATTAAACTCCTACACAGAAAAAACAAGCCCCTATGAATGCGGATTTGACCCTATAGGATCAGCCCGCCTACCTTTCTCTATAAAATTCTTCCTAGTAGCCATCACATTCCTTCTTTTCGACCTAGAGATCGCCCTCCTACTTCCTCTCCCATGGGCAACCCAAACAACAAACCTAAAAACCATACTCATTATAGCCCTCACCCTAATCTCACTCTTAGCAATCAGCCTAGCCTACGAATGAACTCAAAAGGGATTAGAATGAACCGAATATGGTATTTAGTTTAAAACAAAACAAATGATTTCGACTCATTAAATTATGAACCAACTCATAAATACCAAGTGTCTTTAGTATATATAAATATTATCATGGCCTTTACAACATCCCTCGTAGGACTGTTAATATATCGATCCCACTTAATATCCTCACTCCTATGTCTAGAAGGAATGATATTATCACTATTTATCATAGCAACTCTCATCATCCTAAATGCACACTTCACCCTAGCCAGCATAATGCCAATTATTCTACTAGTTTTCGCAGCATGTGAAGCAGCCCTAGGACTATCGCTACTAGTAATGGTATCAAACACATACGGTACCGACTACGTACAAAACCTAAACCTTCTCCAATGTTAAAATATATTATCCCAACTATCATACTAATACCTCTGACCTGAATATCAAAAAATAGCATAATCTGAACCAACACTACAGCCCACAGTCTGTTAATCAGCTTCACAAGCCTACTCCTCCTAAACCAATTCAACGACAATAGCCTAAACTTCTCACCAATGTTCTTCTCTGACCCCCTATCTACCCCCCTCCTAATCCTAACAATATGGCTCCTACCCCTAATACTAATAGCCAGCCAATCACACCTACTTAAAGAACCCCCAACCCGAAAAAAACTGTTTATCACAATACTAGTTACGCTACAAACATTCCTAATCATAACATTCTCAGCCATAGAACTAATCCTGTTCTATATCCTATTTGAAGCCACACTCATCCCAACCCTCATCATTATCACCCGATGAGGTAACCAAACAGAGCGCCTTAACGCAGGCCTTTACTTTCTATTCTATACCCTAATAGGATCTCTCCCCCTCCTAGTAGCACTGATTTATATCCAAAATATCACAGGATCCCTAAACTTCCTAATACTCCAATATTGAACCCAAGCCGTATCCAACTCTTGATCCAACGTCTTCTTATGACTAGCATGTATAATAGCTTTCATGGTAAAAATACCCCTCTACGGCCTCCACCTATGACTACCTAAAGCACATGTAGAAGCACCCATCGCCGGCTCAATAGTCCTAGCCGCCATTCTACTAAAACTAGGAGGGTATGGCATACTACGTATCACAACTATCCTAAACCCCCTAACAGAAATAATGGCATATCCGTTCATCATACTCTCCCTATGAGGGATAATCATGACTAGCTCCATCTGCCTGCGCCAAACAGACCTGAAATCACTCATCGCATACTCTTCCGTAAGCCACATAGCACTCGTCATTGTAGCAATCCTCATCCAGACCCCATGAAGCTACATAGGAGCAACAGCCCTGATAATCGCCCATGGCCTCACATCATCCATACTGTTCTGCCTAGCAAATTCAAACTACGAGCGAATTCACAGCCGAACAATAATCTTGGCTCGAGGACTACAAACACTCCTTCCACTAATAGCCGCCTGATGACTACTAGCAAGTCTGACAAACCTGGCTCTACCACCCTCTATCAACCTCGTCGGAGAACTACTAGTAATCATGTCCTCTTTCTCATGATCAAATATTACCATTATCCTAATAGGAACCAACATCATAATCACCGCCCTATACTCACTATACATACTAACCACCACTCAACGCGGCAAGTACACCCATCACATCAACAACATTACACCCTCATTTACACGAGAAAACGCCCTAATAGCACTCCACATCCTACCCCTCCTACTACTATCCCTAAATCCCAAAATCATCCTAGGCCCCCTCTACTGTAAGCATAGTTTAAGAAAAACACTAGATTGTGAATCTACCAATAGAAGCTCTAACCCTTCTTACTTACCGAGAAAGCATGCAAGAACTGCTAATTCATGCCCCCATATCTAACAATATGGCTTTCCCAGGCTTTTAAAGGATGGTAGCTATCCATTGGTCTTAGGAACCAAAAAATTGGTGCAACTCCAAATAAAAGCAATAAACCTCTTCTCCTCTACCACACTAACAATACTCTTTGTGCTGACATTACCCATTATAATAACAAATACCAACATCTATAAAAGTGATAAATACCCAACATATGTAAAAAACACAGTCTCATCCGCCTTCCTAATTAGCCTAGTCCCAATAATCGCATTCACAAATACGGGCCAAGAAATAATTATCTCAAACTGACACTGAATTACTATCCAGACTCTCAAACTAACCCTCAGCTTCAAAGCAGATTACTTCTCAATCGTATTTGCACCAGTAGCACTATTTGTCACGTGGTCTATCATGGAATTCTCGATATGGTACATACACTCAGACCCACACATCAATCAGTTCTTTAAGTATCTCCTCCTCTTCCTCATCACAATAATAGTCCTTGTCACAGCTAACAACCTTTTCCAACTATTCATTGGCTGAGAGGGTGTCGGAATCATGTCCTTCCTACTAATCGGGTGATGACACGGACGTACAGATGCAAATACAGCTGCCATCCAAGCAATTCTCTACAACCGCATCGGAGACGTAGGATTCATTATAGCCATAGCATGATTCCTATCAAACCTAAACACATGAGACATACAACAAATCTTTATAATTAACCCAACCCACTCAAACCTACCGCTAATAGGACTAATCCTAGCCGCAACCGGAAAGTCCGCCCAATTCGGCCTTCACCCCTGACTCCCCTCAGCAATGGAAGGCCCTACACCCGTCTCAGCACTACTCCACTCAAGCACAATAGTCGTAGCAGGAGTTTTCCTATTAATCCGATTCTACCCCATGATAGAAAATAACAACCTCATACAAACCATCACAATATGCCTAGGAGCTATCACCACACTGTTCACGGCAATATGTGCACTTACCCAAAATGACATCAAAAAGATCATTGCCTTCTCCACCTCAAGCCAACTAGGCCTGATAATAGTGACAATTGGCATTAACCAACCCCACCTAGCATTCCTACACATCTGCACACACGCATTCTTCAAAGCCATACTATTTATATGCTCCGGATCCATCATCCACAACCTAAATAACGAACAAGACATTCGAAAAATAGGAGGCCTATTCAAAACAATACCCTTCACCACAACAACCCTAATCGTAGGCAGCATAGCCCTCACAGGAGTGCCATTCCTAACAGGGTTCTACTCCAAAGACCTAATTATCGAAGCCGCCAATACATCTTACTCCAACGCCTGAGCCCTATTAATTACACTAGTTGCCACCTCCCTTACAGCCGTCTACAGTACCCGCATCATCTTTTTCGCACTACTAGGACATCCCCGCTTCCCTACATCAACCCTCATCAATGAAAATAACCCACTTCTGCTTAACTCACTCAAACGCCTTATAGCGGGAAGCATTTTCGCAGGGTTTATTCTCTCCCACAACCTTCCCCCTATAACAACACCCCTAATAACTATACCTCCCTACCTCAAAATGACAGCCCTAGCAGTAACCATATTAGGCTTCACTCTGGCATTTGAAATCACCCTCAATACCCAAAACCTAAAACACAAACACCCCACAAATAGCTTCAAGTTCTCTACCCTCCTAGGATATTTCCCCACTATTATACATCGACTACCACCCCACCTAAGCCTAACAGCAAGCCAAAAACTAGCATCCTCCCTTCTAGACTCGGCATGACTAGAAAACATCCTACCAAAGTCCATAGCCCACGCACAACTAAAACTCTCAACACTAGTCTCAAACCAAAAGGGCCTAATAAAAATATATTTTCTATCATTCCTTATCACTATCCCCCTCAGCTTAATCCTATTTAACCCCCACGCGTAACTTCCATAATCACCACAACACCAACAAACAAGGACCAACCCGTTACAACAACCAATCAAACACCATAGCTATATAAAGCTGCAACACCCATAGCCTCCCCACTAAAAACCCCAAAATCACTCATACCATAAACAACCCAATCACCCAGACCACTAAACTTGAACACAAGCTCCACCTCTCCTTCCTTCAAGACATATAAAACCGCTAAAAACTCTATCATCAACCCTAAAAGAAACGCCCCCAACACAACCTTATTAGAGACTCAAACCTCGGGATACTGCTCAGTAGCCATTGCAGTCGTATAACCAAATACCACCAACATCCCACCCAAATAAACCAAAAACACCATTAAACCCAAAAAAGACCCGCCAAAATTCAACACCATACCACAACCAGCCCCACCACTAACAATCAAACCCAGCCCACCATAAATCGGCGAAGGCTTTGAAGAAAACCCAATAAAACCAATCACAAAGACAATACTCAAAATAAATACAACATACGTTATCATTATTCCCATATGGACTCTAACCATAACCAACGGCATGAAAAACCATCGTTGTAATTCAACTATAAGAACACTAATGACAAACATCCGAAAATCTCACCCCTTAATAAAAATTATCAACGATGCATTCGTTGACCTCCCAGCTCCATCAAACATCTCATCGTGATGAAACTTCGGCTCCCTACTTGGCGTCTGCCTAATCCTACAAATTCTAACAGGCCTATTCCTGGCCATACACTACACACCAGATACACTCACCGCATTCTCATCGGTAACCCACATCTGCCGTGATGTAAACTACGGGTGAGTCATCCGCTACATACACGCAAACGGCGCATCCATCTTCTTCATCTGCCTCTTTACTCACGTAGGACGCGGCCTATACTATGGCTCCTACACATTCCTAGAAACCTGAAACATCGGAGTTATCTTACTACTCACAACCATAGCTACCGCGTTTATAGGCTACGTACTGCCATGAGGACAAATGTCATTCTGAGGGGCAACAGTCATTACCAACTTACTGTCAGCTATCCCCTATATTGGAACAGACCTAGTAGAATGAATCTGAGGAGGCTTTTCCGTAGACAAAGCCACCCTTACACGATTCTTTGCCTTCCACTTTATTCTTCCATTCGTTATCACAGCACTAGCCATCGTCCATCTACTATTCCTCCATGAAACAGGATCCAACAACCCAACAGGAATCCCCTCAAACGCAGACAAAATCCCATTCCACCCCTATTACACAATCAAGGACATCCTAGGTATCCTACTCCTAATAACAACACTACTCACACTAACCTTATTTGCCCCAGACCTCCTAGGGGACCCAGACAACTACACCCCCGCAAACCCCCTTAGCACACCACCACACATTAAACCAGAATGATATTTCCTGTTCGCGTACGCGATTCTCCGATCAATCCCCAACAAACTAGGAGGCGTCCTAGCCCTAGCTCTCTCAATCCTAATCCTGGCCCTAATCCCAATACTACACACATCCAAACAACGAAGCCTAATATTTCGACCCCTCAGCCAATGCCTGTTTTGAGCACTAATCGCCGACCTACTAACACTCACATGAATTGGAGGACAACCCGTCGAACACCCCTTCATCATCATCGGACAAGTCGCCTCAATCCTATATTTCCTCTTAATCTTAGTACTAATGCCCGTAGCAGGCATTATCGAAAACAAACTCCTAAAATGAAGAGTCTCTGTAGTATATGACATTACCCCGGTCTTGTAAGCCGAAAAAGGAAGCGACACACCTCCCTGAGACTCAAGGAAGAAGCTCAAGCTCCACCATCAGCACCCAAAGCTGAAATTCTAGATAAACTATTCCCTGATTTCCTTGTATGTACTACCTACAAGATTATAAAGTACTTACTTAGTACTATAATCTTAAATGTACATACATACATGGCTATGTACGTCGTGCATTATTGCTCTACCACATACAATAGTACATACTATGTATAATCGTACATAGGACATATTATGTATAATCGTGCATTACACTATCTAGTACATGCTTATAAGCATGTACTAGTGAACTGTTAGCACCACATGGTACATGCTAGTTCTTTATAGTACATGGCACATGTACTCAAATCATTTCCAGTCACCAAGCGTATCCCGCCCCCTAGATCACGAGCTTGATCACCAGGCCGCGTGAAACCAGCAACCCGCTCGGCAGGTTTCCCTCTTCTCGCTCCGGGCCCATAGCATGTGGGGGTTTCTAAGAATGAACTTTATCAGGCATCTGGTTCTTACTTCAGGGCCATCTCATCTAGAATCGCTCATTCTTTCCCCTTAAATAAGACATCTCGATGGACTAGTGACTAATCAGCCCATGCTCACGCATAACTGTGATTTCATGCATTTGGTATTTTTATTTTTTGGGGGATGCTTGGACTCAGCTATGACCTCCCGGTCTTAACTTAGTCAATTAACTTGTAGCTGGACTTTAATTGAACCTTATTTATCAGCCTGGCAGTGAACTACAGGGGTTATTCAGTCAATGGTTACAGGACATAGAAAAGTGCGATGTACATTCATTCGCCACCGGACATGCATTCATACATTGTTCTCTAGCATTGGCAAGACATGCGTATTATTCAGTCAATGGTTACAGGACATACCAATATATATCCCCCCCGTGCTCCTTTAAATCTCTAATTCACCTACTTAAATACGCTTCCCCAAAGCAGAGAGCTATCCCCCTAGATTCTACAAACGAATTTGCTAAAAACTAAATATCAAAACCCGCACAAACCACGCAAGGCAGCCCTACAAGTTAAACTT

ntdensity(hippo)
basecount(hippo)

ans = 

    A: 5363
    C: 4700
    G: 2295
    T: 4044

length(hippo)

ans =

       16402

 basecount reports nucleotide counts for a sequence.
 
    basecount(SEQ) counts the number of occurrences of each nucleotide in
    the sequence and returns these numbers in a structure with the fields
    A, C, G, T(U). Ambiguous nucleotide symbols and gaps (-) are ignored by
    default. Other unrecognized characters are also ignored, but a warning
    message is displayed.
 
    basecount(...,'AMBIGUOUS',AMB) specifies the behavior when ambiguous
    nucleotide symbols are present. Options are: 'Ignore' skips ambiguous
    symbols, 'Bundle' counts and bundles them into the Ambiguous field of
    the output structure, 'Prorate' counts and distributes them
    proportionally in the appropriate fields for the four standard
    nucleotide symbols, 'Individual' counts and reports them individually,
    and 'Warn' ignores them and displays a warning message. Default is
    'Ignore'.
 
    basecount(...,'GAPS',true) adds a field to the output structure
    with the gap count. Default is false, it ignores gap symbols.
 
    basecount(...,'CHART',STYLE) creates a chart showing the relative
    proportions of the nucleotides. Valid styles are 'Pie' and 'Bar'.
 
    Example:
 
        S = getgenbank('M10051')
        basecount(S)
 
    See also aacount, baselookup, codoncount, cpgisland, dimercount,
    nmercount, ntdensity, seqstatsdemo, seqviewer.

    Reference page for basecount

dimercount(hippo)

ans = 

    AA: 1723
    AC: 1451
    AG: 859
    AT: 1330
    CA: 1528
    CC: 1450
    CG: 414
    CT: 1308
    GA: 689
    GC: 666
    GG: 449
    GT: 491
    TA: 1423
    TC: 1133
    TG: 572
    TT: 915

seqwordcount('ACGT', hippo)

ans =

     0

seqwordcount(hippo, 'ACGT')

ans =

    27

    
% Q7: Locate potential protein coding genes in the previous sequence by locating ORFs (you may use
% the Matlab command seqshoworfs.m -- note that it can return a structure of start and stop
% positions for each ORF in each reading frame). Choose an appropriate threshold to identify
% significant ORFs, explain all of the choices you make. Paste the key Matlab commands into the
% report. Show all results. Briefly explain / discuss your findings.
seqviewer(hippo)
seqshoworfs(hippo)

ans = 

1x3 struct array with fields:

    Start
    Stop

orfs = seqshoworfs(hippo)

orfs = 

1x3 struct array with fields:

    Start
    Stop

orfs

orfs = 

1x3 struct array with fields:

    Start
    Stop

orfs.Start

ans =

  Columns 1 through 10

         145         772         907        1744        2071        2680        3100        3283        3616        4309

  Columns 11 through 20

        6199        6562        7807        7945        8116        8362        9169        9988       10186       10348

  Columns 21 through 30

       10795       10954       10996       12064       12472       13225       14215       14770       15322       15454


ans =

  Columns 1 through 10

          47        1232        1391        2426        2750        4022        4472        5612        5687        5960

  Columns 11 through 20

        6587        6947        7340        7694        8795        8966        9275        9707        9923       11141

  Columns 21 through 30

       11642       12251       12650       13094       14105       14258       14570       14981       15518       15596

  Columns 31 through 34

       15692       16004       16145       16205


ans =

  Columns 1 through 10

        1362        2067        2577        2901        3393        5250        5316        6321        7032        7650

  Columns 11 through 20

        8085        8625        9873       10395       10770       11103       11661       12078       12234       12597

  Columns 21 through 29

       12786       13431       14130       15264       15504       15672       15825       15951       16119

[orfs n_10] = seqshoworfs(hippo, 'MinumumLength', 10);
Error using seqshoworfs
Too many output arguments.
 
[orfs n_10] = seqorfs(hippo, 'MinumumLength', 10);
Undefined function or variable 'seqorfs'.
 
n_10 = seqshoworfs(hippo, 'MinumumLength', 10);
Error using seqshoworfs (line 99)
Unknown parameter name: MinumumLength.
 
n_10 = seqshoworfs(hippo, 'MinimumLength', 10);
n_10

n_10 = 

1x3 struct array with fields:

    Start
    Stop

[orfs n_10] = seqorfs(hippo, 'MinimumLength', 10);
Undefined function or variable 'seqorfs'.
 
[orfs n_10] = seqshoworfs(hippo, 'MinumumLength', 10);
Error using seqshoworfs
Too many output arguments.
 
[orfs n_10] = seqshoworfs(hippo, 'MinimumLength', 10);
Error using seqshoworfs
Too many output arguments.
 
orfs = seqshoworfs(hippo, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 100);
orfs = seqshoworfs(hippo, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 64);
orfs = seqshoworfs(hippo)

orfs = 

1x3 struct array with fields:

    Start
    Stop

orfs.Start(1)
Expected one output from a curly brace or dot indexing expression, but there were 3 results.

% Q8: Choose one of the identified ORFs and prove via p-value that the ORF is unlikely to be due to
% chance. Specify for this the test statistic, the significance level, and the null hypothesis. . Paste
% the key Matlab commands into the report. Briefly explain / discuss your findings.
 
orfs(1).Start(1)

ans =

   145

orfs = seqshoworfs(hippo)

orfs = 

1x3 struct array with fields:

    Start
    Stop

length(orfs)

ans =

     3

orfs

orfs = 

1x3 struct array with fields:

    Start
    Stop

orfs(1).Stop(1)

ans =

   247

orfs = seqshoworfs(hippo)

orfs = 

1x3 struct array with fields:

    Start
    Stop

orfs = seqshoworfs(hippo, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 64);
orfs = seqshoworfs(hippo, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 64);
orfs = seqshoworfs(hippo)

orfs = 

1x3 struct array with fields:

    Start
    Stop

orfs_valid = seqshoworfs(hippo, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 64);
 orfs_valid(1).Start(1)

ans =

        3100

orfs_valid(1).Stop(1)

ans =

        3322

orfs_valid

orfs_valid = 

1x6 struct array with fields:

    Start
    Stop

orfs_valid = seqshoworfs(hippo, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 64);
orfs_valid(1).Start(2)

ans =

        7945

orfs_valid(1).Stop(2)

ans =

        8623

valid_start = orfs_valid(1).Start(2);
valid_stop = orfs_valid(1).Stop(2)

valid_stop =

        8623

codoncount(hippo(valid_start) : hippo(valid_stop))
Warning: Unknown symbols 'EFIJLOPQ' appear in the sequence. These will be ignored. 
> In codoncount (line 160) 
AAA - 0     AAC - 0     AAG - 0     AAT - 0     
ACA - 0     ACC - 0     ACG - 0     ACT - 0     
AGA - 0     AGC - 0     AGG - 0     AGT - 0     
ATA - 0     ATC - 0     ATG - 0     ATT - 0     
CAA - 0     CAC - 0     CAG - 0     CAT - 0     
CCA - 0     CCC - 0     CCG - 0     CCT - 0     
CGA - 0     CGC - 0     CGG - 0     CGT - 0     
CTA - 0     CTC - 0     CTG - 0     CTT - 0     
GAA - 0     GAC - 0     GAG - 0     GAT - 0     
GCA - 0     GCC - 0     GCG - 0     GCT - 0     
GGA - 0     GGC - 0     GGG - 0     GGT - 0     
GTA - 0     GTC - 0     GTG - 0     GTT - 0     
TAA - 0     TAC - 0     TAG - 0     TAT - 0     
TCA - 0     TCC - 0     TCG - 0     TCT - 0     
TGA - 0     TGC - 0     TGG - 0     TGT - 0     
TTA - 0     TTC - 0     TTG - 0     TTT - 0     

codoncount(hippo)
AAA - 197     AAC - 175     AAG -  75     AAT - 153     
ACA - 159     ACC - 140     ACG -  47     ACT - 131     
AGA -  63     AGC - 121     AGG -  64     AGT -  63     
ATA - 157     ATC - 151     ATG -  67     ATT - 106     
CAA - 170     CAC - 142     CAG -  74     CAT - 153     
CCA - 116     CCC - 140     CCG -  42     CCT - 142     
CGA -  32     CGC -  37     CGG -  30     CGT -  39     
CTA - 192     CTC - 121     CTG -  52     CTT -  96     
GAA -  73     GAC -  53     GAG -  42     GAT -  48     
GCA -  80     GCC -  80     GCG -  21     GCT -  52     
GGA -  44     GGC -  47     GGG -  28     GGT -  23     
GTA -  52     GTC -  36     GTG -  29     GTT -  42     
TAA - 113     TAC - 111     TAG -  82     TAT - 106     
TCA - 114     TCC - 114     TCG -  38     TCT -  94     
TGA -  72     TGC -  44     TGG -  30     TGT -  46     
TTA -  93     TTC - 100     TTG -  41     TTT -  72     

codoncount(orfs)
Error using codoncount (line 63)
The input structure contains multiple sequences. Please specify a single sequence.
 
help codoncount
 codoncount report codon counts for a sequence.
 
    codoncount(SEQ) counts the number of occurrences of each codon in the
    sequence and displays a formatted table of the result. Codons with
    ambiguous nucleotide symbols are not counted by default. Gaps (-) are
    removed from the input sequence. Codons with other unrecognized
    characters are not counted, but a warning message is displayed.
 
    CODONS = codoncount(SEQ) returns these codon counts in a structure with
    the fields AAA, AAC, AAG, ..., TTG, TTT. 
 
    [CODONS, CARRAY] = codoncount(SEQ) returns a 4x4x4 array of the raw
    count data for each codon. The three dimensions correspond to the three
    positions in the codon. The index in each dimension corresponds to
    nucleotides in the 'ACGT' order. For example the (2,3,4) element of the
    array gives the number of 'CGT' codons in the input sequence.
 
    codoncount(...,'FRAME',F) returns the codon count for reading frame
    F, where F is 1, 2, or 3. Default is 1.
 
    codoncount(...,'REVERSE',true) returns the codon count for the
    reverse complement of SEQ.
 
    codoncount(...,'AMBIGUOUS',AMB) specifies the behavior when ambiguous
    nucleotide symbols are present in a codon. Options are: 'Ignore' skips
    codons with ambiguous symbols, 'Bundle' counts and bundles them into
    the Ambiguous field of the output structure, 'Prorate' counts and
    prorates them into the other codons with standard nucleotide symbols,
    and 'Warn' ignores them and display a warning message. Default is
    'Ignore'.
 
    codoncount(...,'FIGURE',true) creates a figure showing a heat map
    of the codon counts.
 
    codoncount(...,'GENETICCODE',CODE) overlays a grid on the figure
    grouping the synonymous codons according with the genetic code CODE.
    Default is 'Standard' or 1. Set CODE to 'None' to create a heat map
    without showing the grid.
 
    Examples:
 
        codons = codoncount('AAACGTTA')
 
        r2codons = codoncount('AAACGTTA','Frame',2,'Reverse',true)
 
    See also aacount, basecount, baselookup, codonbias, dimercount,
    nmercount, ntdensity, seqrcomplement, seqshoworfs, seqstatsdemo,
    seqwordcount.

    Reference page for codoncount

hippo_random = hippo(randperm(length(hippo))
 hippo_random = hippo(randperm(length(hippo))
                                             ?
Error: Expression or statement is incorrect--possibly unbalanced (, {, or [.
 
Did you mean:
hippo_random = hippo(randperm(length(hippo)))

hippo_random =

CGGTAACCGAAATGCAAACCTAGACCCAACCTCACAAATTCCGAGAATCGATACTGTTACCTATTCATCCAGTTTTCTCCGCACTTTCGACTGTATAATCACACCAGCGGTACCCTACCGGACGCTCCTCGACAAATTATCTACTACAAGACTATTAATCCCACCGATTCATCCTTCTAAAAAACTGTGCTCTCTCGCCCAATCTAAAAACCTTACCTAACTTTTGAATAATATTCGACAGGTCGATAAAAAAAATGCTGATTGACTATTCTATACAATAGGACATCTCCCAGCCGCGAGTACTATCACTGATATCACCAAGTATCTCCTACTACTAAAACAATACTATACATCGGCTCACTAGAACAGAATAATTTAAAAAAATACTCAAGATGTACCGGCTCATCCGACTTGTATCCCAACCTTTTTACCATTTACTAATGTCAATGGTGGACGAACCTCATTGACGAATGCATCTAGCGGACTCATCTCCAATAGCGCCACAAACTCCATATGGTCCCTTCTCACACATTACACCACAGTAAACGATTCCCATTAGTAAACAGCCGATAGTGCTAATTCACCTTCTACTTCCCCTACATAATTATTATCTCCCAACCGGTACCTAATTGCTTCTGGCACATGAGTTCTCGATTACATCACATAATTGTCCCTATATGCGTGTCAAGATCATGCCGGAGATTAGTATGGGGATCGATTAACATCAATCAAACCTCAGGAATCAAAACTCCCTACACTAACTTAGAGACCTTAGCTTTTTAGTGACACAATGCAAAACAACTAATATACCGAATATCCGACCCCACTACCGCTACAATCTAGACATTAATGCACATACTAGTATCCTACAAATGAACAGCAAACCCATTGCTACAACGCTCACAAGAACTGTTTCTATAATGTCCGCCTTCCATCTTTCGTGAGGGAACACCATATCTGTAGAAAGATATAAACGCAAAATAGTCGTTTTACCCAGAGATACTCACGACTACGGAGGCCAACGCAACCCACGCCTTCAGCGATTGATGCCACTATACGGGTTGCATCGGTCACTTAACTTCACACTTCTAAAGTTCCGCCGCAACGCCTCTTCTCCATTGACTAAAGGTCCCCCACCGCAACTCAATAACGTACACCTATCCTATATCCGTAGAACAAACAGATCGTCCCTCTTTCATTACCAACCCCAATTGAACTAACAACATATAGCCCTAGGCCCACTATATCGCTGCAATACGTCCCTAAGCAGGCTAGTTCCAATAATTCAATCTAACTTAATGGCGGTAACACAAAGCACAATTCGAAAGTCAACAGAGCCCCGCAAACACTAATAATAGAACAGTATTTAATCACATATCAAGACCCCCCTAGCTCACTGCCGCACTAAAATAACGGATGATAATAATTGCCCCAGAGGCGCTTACTCCTAAACTAAGCCGTAAGAAAAAGTACTTAACTCCTGTAAGTCAAATCACGGCCGCAAATACATTCCCTCCCCCCTAACCGGTCTTGTCGCCTTTATAACGTAACTGCAAAAATCATCGCACCGGCTCCCTGTCTGGTCTCTTGAGCACGTTACCTCACCAGGACCTGATTTCAATTCACTGACTGCCTATCCGCCTCTTCTCCAAGTAGTATCTCCTAAAAAAATAGAGACGGTCCCCACCGCACTACCGTCCACTGACCTCACTGCCTAAATCCTAAGCACCAAAAAACTCGTTTACCAAGCTATGCTAGCTACCCTTTGTTGTACGAAAACCTATCTGTCGATGTCTAAACGACTGATTAAAAACAACATAAACACTCATCACTTTCGTAGATACAAAACGGATCTACCTAATCAAAGATATACGGCATATGCGGGTATGCGCTTGTCATGTAACAACTGCCAGCGCTTCAAGCGTCTGCGCACAGAGGGCTTCTCATAATTACGAAGCAAGTTACCTATGACCCTACTAGGTAATCGAAGGGACCACATTCTTCCCCACCGCTAAAGAAAACATGTAATAACCGAATCATCAAACACAAACAGACCTACCAATCCGTAACCCTCAAGACTACATAACCAACCGCTTAAAACAAATTAGTATATAAGCTAGAAAGGTGACCTGCATAAACTCGTCATCTCAGCAATCCCGAATAAAAACCGTCTACCATTGTACAAGCTAGCAGTCAGATTCGAGACAAGCATCCTCAAACTGTATCATAATTACAGATTAGCGCACAAGAGCACAAACGAAACTCGTCTCCTAAGCGGTGCCTTAAAGAATCAGCTAAATGACCAGCTACCAAACCCACACCTCCTCCATGCAAATAAGTCATAAACTATAGCCTGTATCGTGATGAACGTACTTACCGGTACAAACGATTTCGTACATACGTAACAGTAAATTCCATTCCCTTAATACCTAAACCATTCTCAATACAACCTCTTGCGGCAAACCTAAAGTGTCCACCTCCATTCGCCTTGTCTCACCAACCTTAATGGTTCCAATCAACAATCCTAACCTCCACTAAAGTAAATCTCAATCTTCCAGATCTACGCGATATGATAGGCAACGAATCTTCACACATTGGACATTCGCGCACTCTTATATTGCCTACCGACGTATTTAAAGTCCATGGATATCCACTTAAAACACCAATCCCTGCTCCCTAATCAGGCGCTAGAAAATAACCTCATAGACCTAGACACTCCTATCACAAATACGCATGAACGTTATTACCTGACCGCATTCCTGACCGGAAATCCATGATTTATACGTTCAAATCGTTCTCACTTGCCCAGTTACACTCCGTACTAGCTACGAGCCAAGAACTAAAACTCGGCAACCGACTCCTATAACTCCCTCTTCTTCTAGAGAACTCAGGTACGCATGGCTAAAAAACAATGCAATATCATTAAAAATACTCTTTAAGCCTGACAAAGAGAGACACAATCCAGGTGACTCGTTACTCCATGTGAGAACTATTAACACTTACCAAGGCTTAAAGTACTATCCACCTATTCCTCAGGAGCCGACCTAAGGTCCACCTGTATTAAATCACCCATTTGGGAACCCACCGGGAACTCGCAGTTTCAACATACCGCAAACAAAGATTCATAAAACAACCAATAACCTCCAGGAATCAATCCAATCACTGTACAAACAGTCCACACGAGGCGAGCAACTATTGATGGTATGACACTACATTTCGAGCATCGGGGCAACAGGCGAAAGTCTCATCATAAGAACACTTACTTTGCCACCACACATTGGAAGTAATTCTAAGACTCCTAGTAAAGGTAGGATCTGTTGTCAACCATTGACGTGACGTGAGTCAAAGTTCAACAACACCCACCCGTAACACAATCTCCTGTAAGCCTACCTGTCAGCTGGCAATCGATTCATCCCTAAATAGAAACGCTAGGAACTTCTCATCAGCTATAACACTTTCGAATGCACACCGGGTGACCTATATTGGTTAAGTACTCATTTCGCCACAACGCCTATAATCCGTAGTAATCTACGTAACTAGTCATCCATATACCGATACTGTTGTAGACACCCGACTCTCCCTCTCTCAGCCCTCCAGTCATGAAACTAAACTTTACGAGAAGCCGACTCAATTGTTCTTCTCCTACTCATGTATAATATTTCTAGATCACTCAGTCCCTTACAAGGTCCCCTTAAGAACAGAGCCATTCACACCTCAGCGCAAATCTTATGAGCCGATGGCATTTCACGTGCAGAGAACCCGTTACTGCAAATCCCGCCCAACTCCGCGGTAAGCATAAACCTCCAATCGTCCTCTCACAGCCCCTCCAAATCACTCCAAGCACGTTAAGCCAGTAACGACTTGTAACACTTCATCTCTTACATCGAACCATAAGAAATTCACATACCCTTGTCAATTACCGGAAAGACACGCACTTGATAGTTATATAGACTCTTAGACGAATACATTCATCACACATTTTTATTCCCTCGCAGATCAACCACTCACCCACACTCCCGATTAAAAATTCCATACCATCCAGAGACCGATCAACCTGACGTGATCCAATCAACCACACCCGTCCCTGACATCCCTATGTAGCACCGACCGGTCAACAATAGGATGTTGCTTATGGTTAAGGACACACTTGAATATCTGCCTCACAGAGGTATACTACTTGTCCATCGGATAATAGATGTATACGAATAAGAACTATCACCACTGCGACCAACACACAACACATCGATAGAAAATAAACACCAGACATACAAGCGGCGACAAATTAAGAATAGCCGACCCCGATAAAAACGGCGAAACCAAAATGCTTCTTCTCTGATTTATTTCGGCATCACATGACCCCTAAGCTCTAAACTCGCTATTAAAGTACCTATGAACTCCGCGAAAATGCTCCCGCCCCGCCCAAACAACCTCTTTCCTGCACCGGAGCCCCCCAATCTAGAAGTTACCGCTGGGTAGCGACTTTAGGGAACCCACCCACTTACCACACACATGCCATTCGTTGCCCACCGTCGTGTTTCCCAACCACCAATCTTCGTAAAGTGGAATATCGAAAGCGGCACTCTATTATTCCTCTGTTTCCGCATAAGCAGGAGGAAACCACGATAAACATAATAAATTTGTTGCACACCACCCACACAATACCCTTCTCACCCGACATTATATGCCAGTAACACCGCCCTAACAGACCCGGGCTAGAAAATACTTTCGTAGAAGCAGACGGCAGGACCTACACGGCCATCGCACCCTGAAACCCCCCAGCTGCTTTAGTACCTATATTCTTTGCGTCAACAACCAAACTATACGTACCTTCATTGCAGAATTTCCGACAGTGCAACCTGAGCAGTAAGAATCTTGTAACTCTACTGAAAACCCAAACGGCCATAAGCCGCTTCCCACAAGATTTTTTACTTCCAAAATTCGTCATCCCCTTACAAATTCCGCCGCACCTAACACACATAATACATCAAATCGAACCGTTAGCGACAAAAGTTGATCGACAGTAGAAATAAAAAGTACTTGCTCAACATAGATTCCACCCAGATAGCATCACAATCCAGTTCTGAAAAACCTTCACACCTACAGCTAAAAACTGGTGGTAGCTGTTGGATAACACCAAAAAATTTGACAGGTTTACCCGCGTAATACAGACAAGACATTTGGTAAACAGACTGACTAATGCTTTTCCAAACATAACAAACTAAGCTTGATGTAACATCTTTCTCAACCTATCCTTTCTAGCTTCTTCCAAGTAACAAATTTTGACCGATGCTATCTCCTCTTAACCAGGACCGAAACTCTTAGGAATATAAGCTGAGGAAAATAACAATCTTACAACGCATCTGCACCAAGGACCGAGCCCAGCGAGAAAATCGTCTTCAACACGCACCCCAATATCATAGACTTTATAGCAAGAGACCCTGCGCCCTCCCTTATTTATAACTTTACGATATCATAAAATCAAGTGTTATAAAACAGATACCATTTTTACCTCAATGTTTCAGACTCACTACACCACTAACGTTGCAATAGAAGATCTTGAAAAAACCTTGCTTCAACTGGATTATTATTCAAAATCAACATAAGCCCATTGGAATCGCCTACGAAGACATATCTCTATCCACCTCCGCCGAATGGAAGCCAACCTAAGTATCGCCACCTCATACCCCTAAAAATGTTAAACATTCCACCGAGGTTGCCCGCTCCATTCTAAGACTTTTAAGTTTAGGCCCATGGCAAGGTATGCTCGCCGCGCCACAACGAACTTGGACTATAACGGTCACCGACGTATATCAAGTCAGACCCAACGACTTCGCGCCCCGTATTGACAACGGAGGCAAATTCCACGAAGCTCAGCTGAGCACACAACAGTAAACGCATACGCTTTATTTCACCTACAAATTCCTGGTCACTCAGGATATCTAACCACATTAATGCCAAAAATACAACCATACTCCGACACGTGAAACTATCCACAACCCCCTCCTTTCACTAACCGCCCTTTAATACCGCTTTAAAGGCTACAATTAGACAATCTACAGCATCATTCAAATGAGCATTCTCCCTGTAATACTTTCGAGCTCAACTACCCCGTCGATCCGGAGTGGCAAATCCTTCCTTGCACTAATACCATTCAAATCAGATAACACAACTCCTATAGAACAACACTAAGACCTTTAACTTAAGCACTCATAACCTTAACCTTCACTGCTCCAAACTTCAAATCTCTTCCAGGACAGGAACTATTGTAAACGAGAGAGTTTAAAAATCCTATTATCTTAAAGCACGACATTATCCTGCGCACTGGAAGTATCCAACTTGGATCCCAAATTAAGATTCCATTCATGAAGAATACTATAAACCAACCACGAAGAACAATATAGGCATCTAAAACTCCCGGTCAACATAGCCGAAACGTCAAGGCCTAACAGACCTTCTCCGGAATTCTCTACTATACCAAAATCCCCTGATCTAAAAAATTAAACTGAGCTTAACAAAATGGCAATGTCTTAGACTTACCAACATCATACAACTTGACCAGAATACTCTTTTCCACAAAACCGATGATAATCATCGATCAGTCCACCCGTTAGACAAACTGAAACTCAATAAGAAGCAGACCGCACCCTACTAATGGCTAAATTATATACTAATACTTGATAACTCATCTTTTAACATCGAATCGACCCGTTGCAGCCCCAAGAATATCACTATCGCCAGTAGCTTCGCCAATTTCTTAACGTTGTAAAGTCCGAAATTAGTATTTAATGGCATCAAGGTTCGAAACACACTCGCTATCCGCGACAAAGTAAACTACTCACTAAATCCTCCAATTATGTCTTTCGACAAACCAGAATTTACGGTGAATCACTCAAGCCCTCCTCACCAACGCTAGATAAAATCGGCATCTAGACGGAGTCTCCGTTTAATAAAATAACGGAATTCGACGTTCCCCACACATACCAAGGGCAACGCTTTACTCCCAAACATTGAAGTAAAAATGGGGCTTCTTAGCATCCATAGCCAGGGCGAGTCCCCTAATACTCCTAGTGAGAAAGCGAAAGGTAAAGCGCGTCAATCTAGCGAGAATACCCAAACTGCTTAGTCACGATCTTTGTCCTACATAATTGACTCCTATGGCATCTTAGGAAAGCATCCGACAGAGGATCTATCCTGCTGATGCATTGCGTTCGTATTAACCAACAATTTACCTGTGAGCTCTAATCTTATCATATGCCCCAACAAAACTTAATCGGTTGCTACATCCAGCAATAATTTTCACCTCACATTCAACTGTCCACGTATGCATTAGCTCTATCCCATTACCCTAAAACCTATCGCCAATCTTAAGCGTCAGACTCACGTAACCTCAAAACTGACTACCTCATTCCAAATCCAAGACTACGATCCTCCCACACATGTTAGCCACGACGATAATAAGCGAACCAACCGTGGAAAATCACATTCGACTTTGTGCTCACCTCTACTGACCATCGACCACGTACACTTCTGGTCGTCAACTATGCGGGCCCTCCGCCAGACTTACAAGCAACTTCTTCGCAGAAACTGACATTAAGAAAAATTAAATACAGGCGCCATACGCCCCAAGTTATCCATTCTTAGATTCACTTTATTCCATCAAACTATGTATTCGCAACTCACAAGGCTCCCAATGCTGACTGACTACTACCCAGACCTTCATATAAACTACGTGCAACTGTAAAATGCTTCACGACCATTCGCTTTCCCCACAGTCCAGAACCTGACCCTAAGAAGGTAATATAACGCCAAACCATCGCTTCAACAAAAACATAATAGCAATTCCAAAATGAACAGGATAGATAAGTTAAAGGTTTACAGTCTTCCCCATAACTATGTTCCTGATCCTTAAATATAATGGGTTAGTACAACGCAAACGCTCAAATCGGATCATTCTTGGATAATAACAGTCCCATATACTTTTCAGAAATACACACGTGAGGGAAAACTAGAATGCCGTAATATAGACACTTTCTGTCGAGAAGAGGTGATTTCGTTCATCCATCCTTAAAGTTACCGCAACCGCTTCTCTAGCCGATTCACAGACCCTTCCCCCCAAGCTTACTACTCCCTGATCACTACCATTAGCTAGTGTACTGTCAGGCCGGATTCATATAGTCACCACCACGCTATGGTCACCCCTAACACTAACCTCCTTAAAATACCAGCAATAGGCCATACTCCACGGATAGATTACCTAAAAATCACACCTAACGCGAACAAACTAGCCTATTATTCTCACAATTGATTATTACACCAGACATCCTAAATTTCTTCTTGCTAATCCCTATCTCTCGACGCTAGCAAAGAGACTACTCCTAACGTTGTTTCTAAACACCTCGCCAATCTCTTAAACCGAAGTTTGACAATACCGTTTCTCAGAACCTCTTCTAAATATGGCAAGCTAACTTTCATCTTAGAGTTATCCAGGTAGTGAAAAGTATACGTTACGTGCTATATCAGTATCCGCAACTATATACTCAAGGATAATCTTCCAGCACACACTTATATACACCTTGATTTATGTATCATAAGAAATACGCCACAGTCTTGAACTGAATCGTGTGAAAATTCTACTCCATGAATCCAGTGACCTTACGAATATCGCGAAATCATCACTGAGCGTTCTTAGACCCTTCTGGCACCTTTATGTTTTTAAAATACTAAACATGCAAACAACTTCTTATTGTAAAGCCTAGCACTTCAGCCACGGATCCTCCGAAAACACCAACTAGATCTGAGTGACAACATGAACTCCCCTTAATAGCGATTGAGCACTCTCATTTACCGACGCAACTTTTGGATCAATTTATCCCCCCCTGATTCAATGACCTCCGCATCTAGCCACCGCTTCTAACATGTGTTAGTGTACAATAAAAGACCCGTACAACCCAAAAAATCTGCTCACTCCGCTATCAGCCCAAGGACATCAGCGTCATGGAACCCCATATCCAAACCTCAATATAACACACCCCTCGTCTCCGTAATAATTATAGTTACTCGGATTGTAAACATCCCTCGCGTATTCTTCTCAGTACAGTAATAGCTTATCTTCCCGATGTGCACCCTCCTGCGCGGTTTACCTAAACCATACTCATCTTCCACCAATCTGGTACATCAATTCCATAACACCCGATGAAGTTTTAACCGCCTGAGTGTTATTCAGCTTCAGCATGTCAAACTTACCAAATGCCGATTGTGGAGCTCACACATTATTGAACACAACTCTTACTGACTTAATCATTCTACCCTACCGCCATAAAGCCTGTTACGAATAGCGCATTTCTAACATAGAATCACCAACCATCTAACGATGATTCTCACCCCATACCCCATCATCTAAGTGCCACCCTCTAGCTAGAACAGTTCTTCCACCGCAACGTTCCGGGTCGTACACTTAACGTTACAAACGAGCCTAAAATTCAGGGTATACTCTGCCCTTTGAATAACTAGTCATATTCAGCTCTATAAACCAAAACATTATTTTAACTTTCCCTAAATTTACGTCAAATAAACATTCCTTATCGCTTACCTCTTTAGTTACTGCTCGACGTACGCATGCGACTACACTACGGGGTACCCAGAAAAAAGCTCAACCACGCAATAACAAACTGTGCATGAACCGTTAAACTTACTAGAAATAGGTCTCATCCCTGCCACTAGGATCCAGACAATCAGAATAGGACACTTCCCCTAAATCAATTGTAATCTCCAGAGGACTAGCCAAATCATGCCTCGTAAAACTAAATCATAGGCACCTCCCTCGCACCGATACTCTCCCTAAAAAACACTATCCAGATCCCCAACTACATAAAACTAAAATCAGACGCCTAAAACTAATATCTTCTCAACTATAAAACCGATGAACACACACCAATCACTACTTTCCTTGACATAGCCCTCACCAGGAATTTTTACGCAAATGAATATCTCTGTACCACTGACTTAGAACGATTCGACACAAGAACACAAACAAATTCATTTAAGTTTTTTCGCTTTACTTTCGTATACTGTGAAGACTCCAGCCCTGTAGACCTTTCCCAACCATTAAATCTAGTCTTACTAGACAGTATAAGACTATGGCATACGGTCCTATAAATACTAATTTCTGAACTCTATCAATTATACAACTCATATTATGAATGTCCTGACCAGTACCGAGCCAGGCGAGCGATGATCAGAGCGTTAAATATCTCAAAGCCATCCCAATTAACGCCCAAGAACCCTCACTAGAAATAGTGCCCCGAGTTCTCCCCACTATAACATCACCACTATACACCTAGGGATCATACATATAGACAGCGCTCAATGGTACTACTAAAAGTCGAAACGTTTTGAGAGAACAAAACCACATGTAATATCCCGTGGATTCTAAATTTCAGGCCAACAGTGTACACAACACGGAAACTACCACTAATATCTTAATAAAAACCCCAGAATTCAGGTCTAAGAAAACATAGATACCTCATCCTCTCCCTCCGTGCTTAAGGCAATTAGCATCCTCTGACATATAGAATTTAATATTACTAAGAAAACGACAACCCAGGCGATTCTATTGTACACTCGACTATACAAAGTGTTTAGCTGACTCTCGCGCACCCGTATCATTCTCTCATTCCCCTTTCTCAAAAACACGTGCCCTTAACCAGCCATACGACAAATTTTCCAACATTTAATCACAAAATGATCTAAAGACAACCATGCATTCGAAGGATAATCCCGCGCCCTCATTAACCTTAACGCTCTCTTACTGATCCCACGACAACGACACTCTTACTAAACTCAACTACAGCCTCTACGCTCCAAGTTAATCCTTACAAGTTTTACAAAGTACCGAAGGTACAATGACCCAGCCTTTGGTATTCCACTCAGACTCCTCCAAATACATATTTGGAACAAGCATCAGAAACCATTCCCTAAATCTTTCTACAGGGAAAGAACCCTAACAGGGAATCAAACACCTCGAACTTTGTCCTCTCTACCGACGGCACATCGAGCATCCCCGATACCATACAGCACTTGCCGCTCATACCAAAGAATCGCCGAACAATCGGGTCCCTAATTACTATTTAGATCCTCCAACTCTCTCTAACTCCCAATAGTAGCCTTAGTGTAACGACTCAAATAACCATTATATCAATAAAACAGTGCCATAGAAAACTATTATACCGCGTTTCAAAACCTAAAACCCTGGACATCCGCGCCCACAAGAAATTCCAGTCCATCCACCTAACAAGTCTACAGATAATAGCTACACATAACAGGCACGACATCAAAAACAACCATTATCCGTGCATCTCAGATATTCCAGCTATATCAAAAAATGAATTCCCGTTCGCACTAACTAAAACCACCTGTTTAATTATGTACGGGCTCCATTCACTTAACAACAGAGACAGAACTCTCCAGAAGGTCCTTGCTGGAATAACATTTCAAGTCCCTAACCATCCCAACGAAGATTCATCAGCACTCTGGCTGAGAATCTTCGAGACAGCATTTAATCGGCCGTCATTGCACATGCGTTCTATCCATTTGACACAATCACACAACGACAGCCCTGAACTTCCAATTCCATAAACCCTCGTGTGCGTTACCTTGTGATACGAGAAAAGTGAAGTCATCGTCAAGCCGCCCTCGATCCGCCAACATCGGTACCGAAGCACACCCACACCTCACGGGGTCGCGATTCTAACCCATCTCATCCAGCTAACTCTAAAAAGTTTAACAAATTATATAGACTTCCCGCTAGTAACAAGAAACCGTTGGAATTCATAGATTTAATATCCATGCCGATTGTAAATGACTACTCATAACCTCAAACAAACGTTTTTTAGCCCTCGCTAGACCGTTAAAACGCCTATGTCTACTAAATCCAAACATTAACTTAGACTACGCTAAAATTGCATTCCCACACAGTCACCTCCTGAGAAGGATTGCTAAGTTACGGAGAGTACGTGTCCCGTGTATGACTCCGCAACAGAGCCTAGTCACCCTGCTCCCAAATCACCCAAAATGAGACAAAGCATAATACAAATTGCTTTCAGATTGTCAATACCGCACCGAAGAATGCTACTCATAACCTATCTAAACTATTGAAACTTGGTGTTTATGAAATATGGCTCAGCGACCCCTAATGCTAAGGCCAAAATGTCTAAGACTTGTCAATCAAACGGCATGAAACCAATCCCCTCTTCAGAATGACCCAGTAGTTTAGGACAGATCCTTGCCCTTTCACATCATACCCTACCAACTAGAAACAGGAGGCAGAGTTGTATAGCAAACACATAAATTGTCACCCAACCATCCTTCTAACCCACTAGAGGATACTGTGCGGATTCGCAAATCCTCCTCCCCCGGTCTGTTAAACCCCTAACTTAGAACTCTTCCATTACTCAGCTCCTCAAATACGACCACTAAACCCCTTATCACGTCTGCAACTTACCGCAATGCCGACCATCCGGCATATCAAATTCAATAAAGACTTGCCCACACCAAACGCGTTAACCTAACACCCTCAACACCTTGCTAATTTCATTAACGACACAACAGAAGCCTCCATAGCCGCTAACGACCCCCGGGTAAAGTAGCTATCGACGTGGGAATAAACACTTGTACCTCAAAAAAGTCCAAATGCCTAACCATTTTCAGTACTTATACTGACAACAAATAGTACTCCCCCAACACCCCTTAAACTAGTTACTAACTTCACGCATTATAATCAGAAATTCTGGTCAACAGATCTCCATATCCTAATCCGTATGTAAAATGCATACCGAACCATGTTCTTCAGGGCACAGTAACATCAGAATAACCCTATATAACTAACTTTCCGTAAGTATGTGTGGTCATCTAAAACCATAACTTTCAATGCGTGTTACTAAAACACTTTAAATACGCAACGGATCCAGGTTACAAAAACACGGAAAGAACGATTCACTACACAATCACGACCCCTCCATGATCATTTTCCCACCCGCACATTCACGCACCCATATTGCAATAGTACCGTTATCAATGATTTCTCAAAACGACCTTCTTGATAGTCAAATTACAGGTCAACCTATTAAATAAACAAATACAATCCAAACACGCTTCACTCGCTTCAACCCAGCAAAATTACGAACCCCGGTCCAAGAATCTCATTATAAATTTCCACCTCTCCATGTTCAACTATTCTAAGTTTCAATGATGCCCTACACCTAACTTCCAAGACTAGAAACCCCAAACGAACCCCCTGATCCCATTTACCTCAGAATACAATTCGTACTCCGAAGCCGACAAAGACAGCCGTATAGCCCAGCTCGCAGAAATGACATCGGCGAGGGTAAATTGACAAGGGTTCTTACATTTCAACTTGACCCCACACCCGCCAATTCCTCGCAATTCGCGCACTCATATTCGTACACGTACCAACAGGGTTATATCCATTTGTTGAAAACCAGCTACTAAAAATAAACCCCTACACGCAGAGACTTCTATCGCCCACCACCTACATAACTGGAATAACAGGCAAATGCAAAACTCTCCGACTTTCGGCATTCGCTGACCCCCCAAATTACCCATAGAATTCTCACTCTGACGAACACGCCGCGCCACGTCCGCACATAACTTAGACGCAATAGGCAGACTCTAGCTCTTACATACAATAACCATTGTAATTAGATATGGCTCAACTGTTCTCGAAGGACAAATATCGCCGTTACCCTATCGCATACTCTACTCGATACCACTGGATTATGAAAAGTGCGCAACGCAATATGAATTAAACCGACGTATTCTCAAACTCACCCAGCCCACACTAATAGACAAGCGCAAATTTAACCTCTTATTCTCCCTTAAAGCTGGTTAATCCCAGCAATAGACACTCTAGTAAACTGAGAACGACCAAGACCCGTTCTGGAACAGAGATCCACCTTTCATACCGCCATATTCAACGCCCGATCACAATCGAGTAGTTTGCACTCGACGAGGTAGACTATCTTCCACCTCCTTGTTTAACACCTATTATCGATGCCGTACCTAAAGGTTCCATAAACTTCAGAAATCCAACAATGACTCTGTATTCCAAGCGCTGAAAACAAGATAGAGAATCATACCCAGTACATCACTCTAATACAGAGCATACGTCGAACCTATAGCCAGAACTTCGGAAACACACTCCTTCGAACCGTTTACCAACGGACGCGTTAACTAACTTCTTTCTAAAACCTACCAGCCTCTCCTATCCACTAATAACGATGAGTATAGATAACTAACGCGAAGATGGCAAACACACAGCCTGCGGCATTTAACTTGGCACACGACACCTGCAATCTATATATCAAGCTCCAACGTCAGTCACTGTCTCCGCTCTCAGTACTAATCCCACACAAACTATCACTTACGAAATAATCTTTGTCAAATACGTTGTACTGACTGGTCTACCTTAAAGTCTCTAACCTAAGTTTATCCCATCCTTGGCAAACTCGTTTTATACAATTAACCCCCACGAAGATGTGATTACGTGGCTTCGTATAATCTCACCGCAAATTTTCAAAAAAATCCTGTCCTAGGAAACTTATATCGTGAACATGCCCGAACTCTCTTCGAGAACCCCTATTACCAAAATTAAACACTAACAAAAGTTCATATGATTGTATTCAGACATCCCTTCTAAACCAAGTAAAAAACCCACTTCCATCAGGGAAACACATTTAAACTGTGACCTACAAACCTTCATACTTGTCGCTATACAAGAGACTAGCCATACACCAAACATCGGTCGACACCAACCAATCCCTACATGCCTTTTAACCCACCTAAACTGACGACAACCGCACACGTCGTCGACACTTTAACTTGTCAATCGCAATAAATACCGAGCTGCACACCGATATCAAAGAAACAGAAGAGCGACTCGGCGTCGCTGTCCCCCTAAGATTAAACCCTTTCTTCTTTAGCTTAGCTTCTACGACCCCATAATGAATTGTTAACCTACTTTTCATTACAGTTAACCGACAAATACTAAAGCGTAGACGACCCGCGCAAGTCACTCCCTACACAACTCTAATACATAACCTTTTAAGCTACCGTCTCACAAGTTGCAGTACCACCTGGATATCTAGCCGGTCTCCTCCGCACACAGCTCGCTCTAACTCATAGTGAATCAAAGTACCGTTGCCTCGCCGAACTCCACAACGGACGTTAAGTTTAGCTCCCCACACGCTGACCCCACCCGCATGGCGCCAAGCTTCAGCAAAAACGCATTCCACCCCCATGTTCAAAATAGTTAACCATTCCTAAGCATAAACATTCTAATGTGCCAAATGCATATTTCAGAAAACGACCCCCAGAATTAGATACAACAACATCACCTTTGAGCTTCCCATAAAACTAGTCTCGTCACTTAAGGAATAGGATCTTTAACAGATCTCCCACGTCCCCAACCACCCAATAATCTGCTCCAAAGAGCACCCGACCTACCTTAAAAATTCTGGTGACATTATACGTCTAAACGAATGGCCGGACACACGGTGACAGTATGCAGAATTAAATCCACACACTACCCTTCTCAAGCGTATAATAAAGTACCGTTTACACCCCCTTCCACAGAGATGGAGCTCTCCGAGCCATTTATAATTGCAGTAGCTTCATCTTTCCTAATGCTACACAGTATACGTCGACGATTACACCCAATAAGTCAATGGAACCTGCCCTCAACGCCTGCCTATATTCACTCCCTCTAGTCCAAGATGCACACCTAAAAGACCACCCCAATACAGTCACACTACTTAATACCACATACCCACCCATGATGTACGTTGCCCCAATTACCACACTCTTGACAAGTGGTGACAAGGAAGTATGTAGGCCATTAAATTTAAACTATGGTCACTCTACAGTACAATAGCATCCTGTTACCCGAGTGGCATATTACTAATCGTGTAGCCTACCCTGCCACTCTCGATCTGTATCTAACAAACATCCATATCACTTCCGTATACACGGTATCGTCCGCGCCAAACTCCAACTCATAAATTATGCCTAGAACACAATATTACTCTCTTTTTGGCAAATCGCCTCTTAATATCACCATCTCATCTTATAAGTCATTCATCAGATTATAGCTACGTACTCCTCGTATATCTTTCACTCCGCGACCCCAATGTAGTATACTAGCCTACAACAAACTACCTCTCTATCATTATCGCGAGGAATACTGTGACAGTTATAAAAACACATACTGAAACTTGGTCGTACCTTACCTTGGGCCTCAGACAAACAACAGGAACTGAAGTTGATCCACAGATTACCGGTCTCAATGTAGAAAGCCTACTACTGCTCAAATACCTTCAGTCGAAACAAA

orfs_random_valid = seqshoworfs(hippo_random, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 64);
orfs_random_valid = seqshoworfs(hippo_random, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 64);
orfs_valid = seqshoworfs(hippo, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 64);
seqviewer(orfs_valid)
Error using seqviewer>handle_struct (line 186)
This structure contains no field named Sequence

Error in seqviewer (line 122)
            handle_struct(seq(i), names, alpha);
 
seqviewer(hippo)
orfs(2).Start(1)

ans =

    47

orfs = seqshoworfs(hippo)

orfs = 

1x3 struct array with fields:

    Start
    Stop

orfs(1).Stop(1) - orfs(1).Start(1)

ans =

   102

orfs = seqshoworfs(hippo)

orfs = 

1x3 struct array with fields:

    Start
    Stop

orfs.Stop - orfs.Start
Error using  - 
Too many input arguments.
 
orfs(2).Stop(3) - orfs(2).Start(3)

ans =

    30

orfs_valid = seqshoworfs(hippo, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 64);
orfs_valid(2).Stop(3) - orfs_valid(2).Start(3)

ans =

        1560

ans/3

ans =

   520

orfs_valid(2).Stop(3)

ans =

        6899

orfs_valid(2).Start(3)

ans =

        5339

valid_start = orfs_valid(2).Start(3)

valid_start =

        5339

valid_start = orfs_valid(2).Stop(3)

valid_start =

        6899

orfs_random_valid = seqshoworfs(hippo_random, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 64);
orfs_random_valid = seqshoworfs(hippo_random, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 64);
orfs_valid = seqshoworfs(hippo, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3]);
orfs_random_valid = seqshoworfs(hippo_random, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3]);
orfs_valid = seqshoworfs(hippo, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 500);
valid_start = orfs_valid(1).Start(2);
Index exceeds matrix dimensions.
 
orfs_valid = seqshoworfs(hippo, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 64);
valid_start = orfs_valid(1).Start(2);
valid_stop = orfs_valid(1).Stop(2)

valid_stop =

        8623

orfs_random_valid = seqshoworfs(hippo_random, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'MinimumLength', 226);
orfs_valid = seqshoworfs(hippo, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'MinimumLength', 226);
hippo_random = hippo(randperm(length(hippo))
 hippo_random = hippo(randperm(length(hippo))
                                             ?
Error: Expression or statement is incorrect--possibly unbalanced (, {, or [.
 
Did you mean:
hippo_random = hippo(randperm(length(hippo)))

hippo_random =

TACACAATTCCACCGATAAAATCCAATGATTGACACTATCAAGGAGTGTACAAAAGCATCCTACACATTTATCAAAACTACAGCAACTGCCTTATTGTCAAAAGTAACATGTGGGTCAGAATTTACAAAGATATTTCACGCACAATACATGACTCCTAGCTTCAAAGGCAAATAACCCCAAGCCATGCCATTTCCTTCTCGGTCGTACTTCCATAAGCGACCTACTTCCGACCGCCAGAAGGAGTAGCCCACGAGCCAAATAACTCCAACCCAAAACTCGTTAAAAAACTCAACTCAGCGAATAATGAAGTCCAGTAGACTTTAGAACGTGATACCACCATTCAAGTAAACACTACAGTCCCTCCATAACAAGCCCCCCGCAACTATTAATCTCGACTACCATAAATTCCGGTTCTTGGTTCCACTCCCCAGCACTCTGATTTCCAAAGACCGATGATTCCTACAGCTTACTTCCCCGCATCTTATATGCCGCCGTGCCATGCCATCCTCAAGGATCACCGCGAAGCATTATGGCCCCAACATTCAGCATCCTAATGGCCTAAGGCCTGCTCTCACGTTGCCTCATCAATCCTTTGATCTACAGTGAGCATCAGGAAGGATAAAAAACGAGATGATGTTATTCTCAGAACCATTGCAAACATCCTTAAGCCTAGTTGCTCTTAAGCCGAGTCTAAATACACGCAACTACTCGCTTGCTTCTTCGGTCTTTTCAATATTCTGCACAAGCCCATTACATTCTACAAACAGTTCTAAGGAGGCAGAGCATGACCACTTCCATCTAAGCCCATAAGCCACATCATATCCCAGAACTCTTAGTCTCTTCTTCTCCCCCCTCCATAGTAGGGTATATAGCTTACCCCCTTAAGTGTGACTTTTAAGCACATAAGTAAGTACTCTCAATCTGGCAGGCTAGGATACTCCGTACCGCACCTAATCAACTGTAGATATAAGAATACATCAGTACGCAGATTCCCTAAATGCACACCTGATAATTCCGTGAATCATACGTATATTTACTGTGACACCAACGGCCCAATTCACCTCTTGCATGAAACTAATACCCTAAAATAACTCCTTCCATCGCTTATCTTAGGTTTCTCAAATACATGTACGTCAACTGTTGAAAGTCTCAGTAAACCTATTCCACACGATGTGCATTTGACAGTTCATACGAATTATGACAAAACGCATTTTACCGTCAGAGTCGACAAAAATCAAAATCAAACAACTTGACTGGCATACCCTGCAATGAACACCTAAAAATAAACGCTAAAACTAAGCTTCAATACTAAGCCTAAGTAAAACTCCACTATTAAATAACACGATATTTTCAGATTTCAAGTCCTAGACACAATCCTCCTTCCCACTCGTTATATATTTATTACTCTACAACTTATATTCTGGCTCCCATTCCCAGCCACATCCAACAGCTTAATACGACCTAGGAAAGCCTGACGATGGCTCGCGAACATACTAGAGCAAGATGCTTCAACAATATGTAAAACTCCACTCCGCTAAACAAATGTTTCCGACTAGAATCCGGACACAGGAAACTACTCTGCTGGAAACGCTTCAAATGGTCGCGGGTTGTAATGGTGCGTCGCATCGGCTCAAGCCACAACTGCTTTCTATTTAGCTCAACCAGACCTTTAGGGATCACTAATTATAGTTTCTTAACTCGGATGGCTACATATGTGAGACGTTACACATATCTGTGATGAAAGGAAAACCAGACATACCGAAAGTTAACTAACCGATGTAACAATCTATACTAAGGATCGTTTCTAATCACATAGTGTAAAAATGAAAGGTCCCCTGCCTAGACTACATAATCTTAAAAGATACCAGCTACTGCTCAAGTCCTGGCTTGGAGCGCTACAATCTGCGACCTGACGATCCATTTAAACTAAATCAGACACTCATGTGTAAGATCCACACCAGTATATTTTAATTCTACCACTTCTAAATCTGATCCCTGCGGCACCTCAACGTACCATCATACATCTAAGTTCATCCCAATTATTTATAAGATACTACTAAGCATATCCGGAACAACAAATAGCCAGACATAAAACTGCCATTCCAAATAGAAAATTTTTATAAACCCTCCTACATATCTCTCTGCGACACTCGATGACCAGACCATTTCACTTCCTCTACCATTTAACTTCACATTCACCGTGCCATAACTGCACACTCTACACAAAGTGACAGTAAAACAATTTACGGCGTCTCTGTAAGTTGGTGCTCCACCAAACTCAACATTCCCAAAACATGTAGCTCAGATCACATTCAGATCTAAAGTTTCCGTTAACTAGACACCTCTTTCCGAACGTCAGTATGTGCACCCATTAAAACTCATTACCGCGCAGAACTAACATTCATTAGCTCAAAACAGATCCATTGTCCAACACATTGCCATACAACCAAACTTTTATGAGAGTCCCTTCACAGGGCTTGACGGGTTAGCGACCTGAGGAAGACCCGCGAACAACAGCACCTTACACCTACGCCGCAACAATACCACGCAAGCAGAGCCTTATCAATGTAAGGCCTACATAGTCTCCTCTTCTTCCCATGTATTTAATGCCTTAAGATCAACGCATATCCGTCTGAACCTACACGGGTCAGGGTGCGAGCTCTTCCTTTCGCACGCTCTTCACCAACGTCCGCAGAACGTATACTTCACCAAATTTACAGACGTCGAAAAACGACCCGTTAGACAAATCTAGACAGAATAAAATGCGCCCCCTTTATTGTCGGAAGGTCAAATAGTCGCGTTCAGTACTCTATCCTCTCATCAGATTTTTACCACAAGTCCCAACAGGTTGACAGTCGTCGGCTTGACTACTATTAATTAGTGCATATAAACCCACCTAATCCCCCTAGTCCTATTATGGAACCCCCTTTACAAGAAACTTCATCGACTCACTCGCGAAGAGCAGAGATCTTTCACTTAAGATTGATCGCCAATAGCATCCGATAGATGCCATGACCATATAACTATAATCTGACCTATCCCAAATAACCCATACTGAAAGGCAACATTATAGCGAGACAGAACCTCAAAGACACGCCCATAATTGGGAAGGAGATTTGAAAAAAAAAAAGTAAAGCAGTCTCAGGTCCACCTTCGGGGACCTATTCTCCGACATCCTCGGTCATAGCCCGGCGAGCGTTTCTACTCCTTAACACGACTACTTCCCGACTTAATTACGACTTCGCCACAAACCACCGCGCTTTCTTCTCGTTTGGCGCCAACCCCCTGATACGATCGGTGCTAAATCACAATCAACTCTGTCGCTACGGCATACCGACCTCACAACCCCATAGCATTTCCCCTTATCGCTTAACCCCTGCCAACTGCCAAACCAATGTAACCCTTAACAATGACGAAAATGTATCAAACAAAACGTGCCTAAGCCCGTGAACTGCAATCCTAAGCGGATGGTCGTTTCACCATAGTGACACCAGAAACCTCACCATACTCCAAGCCAACCACCACAAATCTGTAGTGCACTGTCATCTGCAAGACCCCATTATGGGAATATTCCTTCTTTTAATTTCAATCAATCAGGTCACATGCACCGTTATCCCCAATGAACGATCCCCCATTACTGGTTCACTCACACATAAAATGACGCCCAAAATCTCCACGATAACTAACTATTTAAATCGACGTTGGTACGTACCAGGTTTAACGCAAACTGCTTTCCGAAACACCTTAATGGCAATATTTATATGACATCGCAATAAAATGGGAGTTATCCATGTGAAACGCTTACTATAAACTCAATATGATCAATCAAATAGCCACTCACCATCGTAAATAAATGATTGTTGACGGAAATGATCGAGTTATCATAGTCCGCCTATTCCCTAACGCCGCCATCTGGATTTCTATCTAACCACAACCAAATGTACCACAATTCGGCTACAAATCCAGGCCCAGGTCAAACCTGAAACCAATTAGCATAGACAATGTATTACCAGAGATACATCCATCACTCTAAACGCACAAGAACAGGTAAGGAATCGTCACAGCTAACCCAAAATCACCCCGATCTAGCCCCTCCAGCAACCACTAAAAATTAAGCTCATCTTACTTGATTCCATCTACCTACACTAAGTATTACCACGTACAAAACCACAACCCACACCCATCAACCAGAATCCGCACGGAACACAGCACAAATTAGCGCGCTCAAAAGATACACTTTGCAAGATACAAAGCCCCCGAAGTTAATTAAACCTAAAAGTTTTTACTTTTATCTATGCGCCGCTTCGAGCATCAAACTTCGCAGCCCAAATCTCCACTCGACCCTCCCCACACGGAATCACGTCTTAGAACCTGAGTAAATCATCGCTCATAACCAATTCGTTATAAAGTGTACTCGACTATAGTCGACACCATGCTTTCGCGTGAAATACACCCTAAAAGGGCAGTCCAGCTCGTCCCCCGTTCTCAACCAAATGACCTTCCCGCATTAAACCCATATCCATCACCACCCACCCACACCATCAAAATACATACATTGTACGTAATCATGCGGACACTAAGCAACCACGAGAGCGCACAGATCAAGTTAACATTTTCTCACAAACAACGCGGCGCATCCACCCCACCTAAAAGTGCTCCCCACCCCATCAATATACGGCCCAAATAATAACAGTCTACACAGAGACCGGAATACACTCCCTCCATATCCTCTCCCACTTCTCATCTGGTCACTAATTGAATAACTTGTATTTACTAACTCCTGCGCAGAGCTCTTATGATACCATACAACAAACAGATATTTACAATAGAGAATTCACCAAACGTACTACCCAGCGCCCTAGGAAATCATAGTGTTCTTACAATAACCACAAATATCTAGACCTGAACAACATATCGAACTACAATGCATTGCACGACAACTTGCGCAGATTAACTCCCCCCCATACGAATCAACCCAGTCAGGATCAACACGACCTCATAAAATTCTCTGTTCCTCCGCAATACGCATCATCTCACACAGCCCTCTTAGTCTCCGCCCGAATCCCAAAACTATGTCTTTAACCAAAAGAGCCGTTGCTCATAATATTAAATGCTTCGAGTTACCAAAACACCATACACCACATTACACTGTAATAATATTGCACCTTGAAAATCAGATTCACTACCAACCGTGTTTTTGCCATAGCAACATATATAAATTGACCAAACTATACTTAAGTACAAAAACCCTCACCGAGAAACGATTTACCGTTCAGCAAGCCACTACGAGCGCAAACTTCCCATATGCGTGTACTACTCACTAACAATGACTTAAACAGCTAACGTATTACAAACAGGATCTTTGACAGTACGGCAAGTTTTAACATTTACCCTCTAATGATGCGTACAAGCTACGCATGTCGAACGTCCCACAATCTATACTAACACATCTATGTTAAGCTCCATCCATTTCCCTATCGGCCAGTATTCAACGTCCTCTACGTCACCGCTCATGCTACGTCTGACTAGCCGTCAAAAACATCACACTTAGATTTCAAACATAAAAGTGCCCCAGAAGCTCGCCTTTCATTGGCCACAGACCTCATCAAAGAACACTATATGAACGCTAGAATAGTCTCCCAGCGAGGTTTAAGCCCCAATTGATATTCCTCAACGTCTTGAATGGCTCAGTTCCATCCCGTTAATAACCTACAACACACTTTCTAAGACCAGAGACCACCTGACTTTAGCCACCAACTCTACAACAATACACATTCTATAGTACGTCTTCGAACCTACCGCCATATTAACAACAAATATTACTCGTACATAATTTGGATTTGCTGTTTAATTAAGCCTCTGAACTCAATGGCTCCCCCCCATTTTACTGCCGTGAATACGGCCGGGCGTTAATACAATTTCCTGTTATCACACCAACAGGTTATCCCTGCCCCTAAGCCCCAAAATACGACCGACACCAATATTACCAACCTGGCTATGAAACTCTATCTAGTAGTACTCAATATCAAACTAAGTGTCACCGGTATAGCCATCCTATCACCCGACTCAAACTACTATAACGCAACCTTTCTCATAATTATAAAGCACAATCTTGTTATAATACCCTTTTTCAAAAAATGCGTCATCTAGAAACCAATCACCAAAGCTAAGCAAGGACTCCCCCCATACATAGTACACTGTCTGTCGAAGACCCGTATCGGATTCTAACTCCACTCAATGCGCTTAACATGCTAAAGTATCAAAGTCCCACGTCCTACTACCAGCAACAGAAACCATCCTTGTCACAGCCAAAAAACAATAGTATTCGTGACTAAGACCGACTAGATAAAATTCTCCCAGGATAACCACCCCGTTTCAACATTACGCAATTACTGAGAAACAATTTCTACTCTAGCATTAGTCCACCATAATCGTAATCACTTGCTTTAACAAAAAGTATCGTCAACTCCTACTGCTTACCTTTCCCTCACTGATCAACGCCGAATCGATTCCCCCTGTGCAGTAACCGGAAAAGCGAATCAGATCATATAGAATCTCAGTACCTAGAACACTGCCTGGTAAACTAAAATAATATATACAAAATGTTTCCCAGAAACAGTAAAAAGAACCATCCAATTACGTAGGGTTACAACAACAACGACACTACACAATTTGAACTCACCAATCTCTTTTTAAAACGCATATTTCAATAGCCGTCACCGTGCAGAAACCTTACTTCGAGCCCTCACCGAGCTTAACAATGAAGATTGTTATACAAACTTCGTAAACGGTTTCTCAACGCACAGTGGCTGAACCTCCTCTTGAACACAGCTTGATTACTTATATACTGCCAAACTCCTAACTCCAACAATTTGTTTTCTCACCTCCCTATGCGATCTAATAAACCCATAAACTTTAATAAGAGTTCACTCAACGCAACCTTCTTACACTAGCTCTCTCCCCCTTCAACAGTCCTCGCTTTATATCAAATAACCCAGACGCCGTGGAGCACATTTGTCCATCAAGCTTTCTTACTTTAACGCATTACACCCATACCAAGAACTCAGAAAATATTACCAACAAATCATCTGTCCATCCTATCACACTATCCTATAACCGGAGCCAAACATATCTACTGAAAGAACCTTGAGCTAACGGCAAAGGATCAGGAACACATCCCACTCTCTAAATACACCGATTAATCCAATCCTGCCTCGATATTCCTGTACACATCGACTAAGCTCCTTTGCACAAAAAACAGAACCTTACAATCCCATTCCAGGAGAGCCCACCTTCCAGATGCATAATCGATGAACTCGATTGGGACTATTACAATAGATGTTCTCCACCCAAGAAAGATTTGTTACCTGAAACACTTTTTCTTCCTTCTGCCTTACCCGAATACAGTTACTACGTTTGTATTGCCGTCTTTATCCAGCTCCTGATTAACCTTAGACCCTGCGCTATAAATTGTCACAAAACCTGATAGGTCTATACACAGCCACCTGAACGCCACGTCATATCTAAATATAAGCATACTACTCACAAACTATGACAATTAGCACATATTACCCAAGACGATTGCACATCCAACGAGAAATCACTACCTCCGCGGGAACAGTTCCGTTAAAACTCCTACCAATACGACTTCTACCCCGCCCTTTTCTGGTCACAGAACCTCTCGCTTAAAAAAATTCCCGCTCTTGACTCACTCCCCCATCTGATTTCCCTACAACCGGATGCCTATTTTCGCCGTTAAAATATCAACAAATCCTAATTTCTCTTTCACATCATTCCCCACAACAAGTGGATAGCGCAGCAAGGCCGTACCTCGTCTCCCCTTTAGTACACTACTACACCCCATCTATTCTAACACAAAGAGGAATCAGCGCAACTCACAAAGACACAAAACGCTGCTGTACATACCAAGCCAGAAACCGTCACACCACCCGTTAACCATTTGACTCACTACTACAAAATAAACCCTGGCACACATAACAACGCGTCATTAGATCTTAAAACGGTACATAAACAACTTGTATAAGCCTAATACTTTCCAACAAGCGATCTCCCAATCAAACTTCAGTAATGCGTCATGAAAAATTCGAAACATCCTTAACAATTCTCCACTATATATCCTCTGAAAGCCAATACATATTCGTTAAAAACTCCGCCACAATACTCATGGATGGTTCAAGCTTTCACGGTAGTTATCATACACATACTCACGTTATTTCGCAGTCTCAATTCGATACGGATATGTTTCCGAATCGGGCCTTTAGCGCCGTTGAAACTCTTTCAACTGTGGATTTTCGTATCAATCGCCCCTAATATAAAGCCCACGATTCGACATGATTTCCATCCAACACCAACATTAGCTGACACCTTGCGAGCGTTAATGTTCACGCGGAATAATAGAATTGATAGTTGCCTGTTCCCAAGGCCCAGCCTATCCATGGACGCAATCCTAGAGTACTTATTTAAACCGAAATTGACAACACCCACAAGCATCATGCCAGCCCTTATGACTATCTTTTCCGCAACCAGAATAGACCCTCCCCAAAACTAACCAGCGACCGAATCCTCAGTATAACTCCATCTCAAATAACAGAGACCAAACAGTCCAACCTTCAGTAGAGGACCCCCGCGAACGAACATATCCCTCATTCCTCGTTTGCAAAGTCAGGAAACAACAATGCACAGGTAACGTTTCCAACCGGGCAATATCCACGATCCGAGTATAAGATCTGGAATACCTTCTGGAAATATTGCCTTGACACAACAGTATATACTGTGTCTTTACGAAAACTTATCGCTAGATTCGCTGAGGTAGCTATCTTATATCTTTTTTTCCGCTATAGCTCACCCCCCAAAAGCCCGTACAATAACACTCTCACTCCCATAACCCGCACACAGCTACTTCACATTTAAAGACCAAATTAATAACTTAATAATATATGTAAACCTTACTCGGTTATGACCAATTCAACAGCAGCGACCCTATTACCTATTCCCGGTTCCTTAGAATTTGACCCATCCGACATAAAACAGAAACGGTACGAACCTTCTTTGAGTATAACCCCCACCGCCCTTACTCAGAAACAACGCCAGGAGCACACATCTTCCAACGCCATCCTACGCTCGTAGCGGCGTTCTGCTGGTCTCAGGCGACCGGGCATCGCGAATCACACAACTATGATTAAGAGACACCAAATTACCCTCAGGATTATAAATTGCGAAGATAGTAAACCTTGCACCATCATTAAATGAATCCGCACAGTTCCAATCCACCTCTCACTTCTGTTATTTGCCACGCCCCATCATTCAGCCACACGGAAACGACCATAGTTAGTCCTATCCGCAGCAAATCTTGTTGAAATCACGACACGTTGTATTAACATATAATATAATACACCACTAACTAAGTTATCCCCCGATGATCCAGAATACGTTAAGTACTACTCCGATCCTGAGTCCGATATGTTCAAGCCACGCTCTTACCGACATTTCCAAACCCGCTAAACTGTTATACAGGCTTCCAGCGCGCCCTTTGCTAGCTACAAGTAACATTACATCTGACTCGTCCGCCATCATAAACTAATACATGATGCAGCATGCATTTATACCTAGTCATCACAGGCACACTGTGCTACCGATATCAAATAAAAGCAGAAGTTACCGTAAACCAACACAATTCCAAATCAATACGTAATCCCTTATGCTAACTAACACCTCAAGTAATCCAACGAACGAGATAGCAAACATCGTATAATCGACTATAAATATACCCGGATGCCTACATCCATTCAGACTCTCTCAATAACTCGACTAAACGGCTGGTCACTCATATCCTTGTACATTCCTGTTTCACCCATGGTCACACCGAGTTAAACAGGCGGACCCGATTACCGATCCCCATCACAAACCTCGATTAACGTAACCGAGTCTGCTAACGAGTACACCAAAGGGATTAATCAACCATAACCCCCTTAACCTGAAACACACCGATTACCGCCCGGCATTAACCCATACCTTTCTAGCGCACCTCACTCACAATTTAGTCACATCCCATTGAGTATCCTCAAATACTTATCTTAGACTTTGAACCTACAAGAAATCCTAGAATTCCATAATACACAGCATATCCCAAAGCCACTAACACGTAATTGAGTTTTCCTCAAGCTCCATTGAATGACACTGAATTATCACACCGACACTAATCTCAATGCTTCACTAAGACGCTCATCTATCAAATTTTACCTACCAAACTTAATACATTCCTGTTCCAAAGACATTGTTACGCCCCCGGCTCTTCAATTCGCACCAACAACTAGCTACACTAATTAAGAAACAAAGCTGGATCGCATAACACGGACAACCCTGCAGTAGGACCACTTAGTTCACACCGCGATTTAACTACAGCATAACACTCTCTCACACGATCGCTAAGTACCCCCTAAAACCCACCTTATTCAGAAACCAACATACCCGCTAACTGTGTTCCTAAAGACTCAGCGTTTCCTGATTCCGTCGTCGAACCTACCAGGACTCCTATTCTCAGTGAACCCTAACCAGCTAATTGCCTGATTGATCACGTTCCAACAACACCACATCATTGAATCTCAGAACGGCCCAAAGAACAGAACAACAATAAAACTCAGAAAGTAGAGCAACGCCGTTCAATCCCCACTCCACAACCAATTACGACAACAGTAATGATTGGTAGCTGCCAATTTTTAGATCTTCTGAAATATAATTCAGTAAATTCCGAGTACCTTACATACAGTGCAATCAATCCTGTCCAACGAGCGCCACCTGTCCACATCAAGTCACTCATTAGCAGTACCTCCAGTTCATATGCAGCCCAAGACATACAATAATTAATCCCATCCTAAGGGTCATCTTAATCACCAGATGTCCTCTCTACAGATCGAAATCTTTGTGTTCGAAAAACTTTAATAATGTAATGAAGGCAGTAAACACGAATCCAGGACTAATTTCGCCACACAACCTATCAAACCAGCATTAGAGTACTCTCACCCGCTTCAGTAAACTCTCCACCCCACTATACCGATTAAAGCGAAATGAATCCGTACCATGATACCCATCTTTGCTCACAGACTCAACTCAACACCACTGTCATAATACGGAAAAATAATAACGTAACATTCCAACATCTGTGCACACTAAAAGCAGAGCTACATGAAACAACTACAGTAGCCATTCATTGGAACAAGCCTAAAAACGTCTAAAGCAAAACAAAGGATTCTAATTATACCGGTAGGAACCATACCTCACTGATGGTTCGTACTTGAGCAATCGCGTTTCCACGAAAATACCCACGTGGAGCGTTTACAAAATGGGCTTCGTTATAACATCATTGTGGTCGACAAATACATCATTCACCCTGCATATAATTGGCCAATGACGGATCTGCCTGACTCTTTAAGGGGACTCTTTCAAAACCTTCTAATGGCGGAGGCCATAAAAAACGTTACCCCCATGTTAATATACCTATCCGCGCCAGTTCCCAAAGATTCCGTTCCACAGACCAAGAAGTGTACATTCATCATATACAAACACCTTTGATTTAAAACAACACCTGGATTAAGGACCTCCCCAGGATACTCCAAGAACTTACACGCATGCCTATTCTTACTCGGGCGAGTCTGAACTGGCTAGTTAAATAAAAGACCAGCCCCCATACTGACCAATAACGGCCTGTACATATCTCCCATCTAACCATTTTCATACCCCCATGGAAAGAATCACCCAAATGACTTACCCCGTCCAGGTTCCTATCACGGCATCATGTAATTGCGCGCCGAGCATACGCATCTTACCACAATAGGGCTTAACCCAATTCAGGTCCTCACATCAACTACAATACTGTAGTGTAACCTACTATATCATTTTTTTCCTTTCTCCTATACATCCTTATTACCAACAAATAGTCGTGACAAGGCACCCGACCAGTCGCTAGACGACAACAGACCAAAAATCAGTCCTATCATCGAACCACGTTAGCTTTTCTAGCGATTACTATTATCAGGTCCACCTAGTGCAGTCACAACTAAGAGTTCACCTACGAAGAATCCCCTAGACCCAACTGCACATACTAAAGCACGTATCTTTGTGGACGGACACATGAAAACTCCAAACAAACTCCAGTCGAGTACAAATGCCTTCCGCAACCGTGACAAATCATCTGGTCAAATATACGAAATATAACTCCATGCTCAAAGTGCTCTCCTGTTGCCCCAATCACGCCTATGGTACACACTGATACATTCATATAAGACATCGTGTCTTAGTTAAAAAATAGGTCTACTGCAAAAACACATTCCCTCAAACTAAAACCTAACTTCATAATTGGCCATGCAAATCACTATAATAGACTTACCCTAAGCTGAAGCATCTGCACACTCTTATACCCGACAATCAGTGAATTTATCAAATACTCGACTCTGCACGCACCTTGATTCACAACAGATAGATGGGTTTATTACAGACATCATATTACTCGAACATCCGTCATCGTGCACGTAAACACCTCTCTGGATCGAAATATTATTTGCGCCCATAATACCTGGGCTAGTCTACTCATCAACGCAAACAAGGACACAACTGTATACCCCCAACGTGACACCATTGAAGTGCAACGTACTCCGTCGTTCCCTGAATTCTAGTGGAGAAACCATACCCCGAAGCGTCAAGTCGCATACAATAACAACTTGTCATAGACTAATTAACTACCTGTATTTCTTCCCTAGATTATATAGTGTCAACAAGCCTAACATTCTACTACACACGGCCAGTCAGTAGACCTCTACCCCACTAGTTCAAAAAGACCACAAACCTGCAGACCTATTCCAGATTCCATCACTGCGTCACCTTTCGACGGTTCTCTCACACACACGATGCGTATCACTCCCCTACCCAGTTCTACTGAATCGAGCAATCTAATCGCACTGTTAGACAATCAATAACACACGTTATTATTAAAAATCGCACAAATTCCACGTTGAAACAGCACAATTCCCGCACTGCACGACACCGTAGACTCATATCCGCACTGAATTATGAACTGAGCACAAGAATCTATAAATTCACGGGATCCCCAGCATGTCACGGTGTTAAAGTGCGAGCAGATAATCCATACTGGCGAGCCCACACCAGTACATTTCTCGCTCACAACAAATACATCGCAGGATCCTCTCCAACCGAGGAAAAATAATGCAAACTCGTATTAAGACTACCCATAGTACGCACTGGTTACAGTATCTCATTCACCAGCACATGGTATTTTCTCGTCGCCCTTAAGCACACCAGATCAATACCCCACGGAGCTGACAATATGAAACCCAAAGCTAATAAGAACCTAAATGCCCTGACGGCATTCAGACGAATGCCGTTGCACCCAGTACTGTCTTCACAATAATTCATAATAGCATCGAGGCTTAGTTTAGGTTATGCCCGTTTCCATGCGCAACAGCCGTAACATAACTCCTATATAAAACACATTCCTTTAAGCCAGCAGATCTTTTACATCCACCATACGCCGATCAAAAGACCTCTATTAAAATATAGCTATAAAAAGTTCATTATCAACGACAGCCCTCTTGGGATGCAAATAGTCCATTCCATTCCATATTTAATTGACGACCCGCGGGTAACGAACACTCCAGCACATCCCAATCAAGTAAACATACTTACCATCACACCTACATCAAAAACACGCTCAAACCTGATGTCCAACACAAACACAGTATACACCCAACAGCCCGCAATCACGCAACCTCATTTAGAAGGTCGGGTCGTTCCATATAAATGCTCAACCACGCACCGCCAATAACACCTCAATCGCACGTCATCTGTGCTCAAGCCTAGGCTCTATTACAAGTCTCTCCAATATTTTTCGCAGCTATAAGCTACATTTCGACCTTTTTGAAGACATCCATACCACTTATACAACAAAGACTATGAACGCTTTATGGAACATCCCCCTCAACGTAAAACAGTGTAATCGCAATACACTTATCACTAGGCCGAGTGTACCGAATTGGATAACTAGCTCTATATGAAACTAAGGTCCTGTATTATCATTCACACTGGATACTTAATCCCTAAGGACCACCCTAGTGCACCGGATAACGTGATGAGCGACATCGCTCTACGTAAGAAAATAGTATTTAAACTTATACCACAAACAACAGCACGATACCAAGACGTTTAAAGAAAGATACTACAGCAATATATCCTACCTCATTTCCACATGTCTAATTTGCAAAACACTCTACAACTAACCTACTATAGTGCTCTTTAAACCAACAAACAATCTAAACTGACAACCTAAGCACGTATCCGTTATGAGAGCCCCATTCACCGGCATAGGTTAAAAAGTGTAATTCACACCACTCTCGGCTTACGCCCCAGTACAAGCCTGCACTCACAACTGTTCACCTTTAACGCCAACCAACCCCCAATAACCCCGCACATAGTTCTCCAGACCAATTATAATCCTTCATCACAAACAACGCGCACTGCAGATTCAAACAACCTACCATCAGAACTTCAAATAGATATTATCTACGCCATGCGTTGCTACATTTTGCGGTGAAACAATATCCTGACTGTAATCATGAGTAAGCGAAGTGTATATCCAGGCATATTATTAGCCACATCCTTAGGCAATCGTTTAGTCTCCCAGTCTGTCACCCACCGGCCTCGAATTCCAACCCGCTCTAAGCCAATATGAATTAACGCGTCCCCATCCAATGCATATACAGTCGTCATCTATTCAATCAGATATCTTCCACCAACAGTTCCTCCTGCACAATCATTTTAATGAATTATTCCGATGTTAACAGTATAATTATTAACTAACAAACAGAGAATTACCAACGAACTGACAGTAACCCCAACCTGATATCACCAGTGCTTGCAGTAGGGGACCAAATAGACCACATCACCTATCAACGCTATGCCGTACCGGTTATTAAGGTTTGAACCACATTACCAGACAGATTCCGGCACACAATAGATCACACAAAAAGCAGTTTCGTTCCTTTAGCAGCAAAACAGAAAATAGAAACCTCGATCCACTTCCTTTGTCCTTCAACTCATAGACCCATAACTACCGATGCGACCAAGTGTGCTAAATAAGGCGACGTATCGCTGCATATCTTCTGATCAGAAATCAAACTCCAGATAGACACATATAATCACGGATACACCCATACTAGAACACAGGACTGCTGACACCTGGAGACTATGGCTATTCCGACCGCAAGCGATTCTGCAGCGTCGTAAAACTCACCTATAGATCATCTTGCTCACCTCTCATACAGCGTAACATTGGATAAATTAAGTTCAGGTTGATATTACCATTCATTAATAACTAAATCCGCCACTTGGCTCTGCGCGCAACTTGATTCTCAGTTAGTATTACATCGCTGGTCCTGTGGTGCCCCGCTCTCTGTAGCCAGTAACCCGATTCTACTTATTTACCCTAAGATAGACTCTCCAACTAGAGTCTCGATCTTATGAAAAGAACCGCGATTTGTCCCCACCCACTTAATCCGCATTGTACGAAATAATCAAACCTAGGTTATCAAGAAACTGAACAATAACCGGGGATTAACCCCCCTACATCACATACATATCCAATCAAGCCAAAAAACATATACTTACTTTGCGCGCTACTACCTCGTTTCCATTCAGGGGCAAAAATCCAATTCCATATTGCCTAGATCCCAGAAAACTCTTGCACCATCATGTCCTCACCAACGCCAATCCTAAGGCCATAAAAAGAGCCTGTTACCACACTCCGTGGTCAGCCCGATACCTCTATTTACGACCATAGATGAAGCGCGATCGATGACGCCAATCACTACCTTTAACCCAGGCCACACACGGCGGGGCACTCTAGGAGGAAGAAACCTATACCCGGCTATTCCTATAACCACAAGAAAAGGCAAACGTGAATGACGCTGTGCTAAAGACACCTATATCAATCGAACAACACGGTCTAACAGCACCGGTCACAGCCTATTCCACAGAAATAATATCTGGTCAGAAAACCATATACTCATCTTGACCTCAAATCCCAAACTACAGAACTCCATCCCACCTGTCTCGAAAATAACGACTTAGTAAGTAAAGAACTGGTCAAATAACGAATCTCGCGCTCGAATAACAACAGCGCATAACATCATCACTCATATAGATACGCTCAGTCACTCAGACAACACATAAAATTTCGATGTTTTTTCTCATTTGGCAAAGATAACTTAAACTCACTCGCACCCCCTATAAAAAAAAATTTCAATACTCAGAACTAGCTCCTCAAAGGTCGTGACGAATCATTATTATAACGTTCTTCTCGTCACTGAATACCTC

orfs_random_valid = seqshoworfs(hippo_random, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'MinimumLength', 226);
hippo_random = hippo(randperm(length(hippo))
 hippo_random = hippo(randperm(length(hippo))
                                             ?
Error: Expression or statement is incorrect--possibly unbalanced (, {, or [.
 
Did you mean:
hippo_random = hippo(randperm(length(hippo)))

hippo_random =

ATATAGGCTACAGCCCAACTATCAATAGGCGAATCCATTCTAATCACATAAATGACAAATGTAATCACATGCACCCTATAAGGATCAACTTTGAACATAATGGCCCTTGTAAGTTCGTACACAAACTAGACTAGGCATGTTCGTATCTAGACTTACCTCTTCGGGATATCAAGAATTCTATATAGCTCTTACACACATCAAATTCACCACCGTTACATATCCTTAGGAGAGATCTATTATGCCGAAGTCCCCCACGCCTTGTCCCGTAACTAGTACGGGCGACTTTAAGTTACAACAGCCACCAACAAAACGTTAAATACAATCCTCCATAATCTAGCTTCGTGTGCTGGTCCACGCCCAAAATAACTCTTTTCGTACACATAGCTAAAATCCCGCAGGCAGACAACTATATGATTCAGCAGACAATTATTCAGCCACAACAAGGTACCTCACTGATGCACAAACAAGAATACAACCCAAGCCCCTTCACATCACACTCCGTTGACTAGAACAGTACATCCTCCCCATATAAGGATTGAATTCGGCATCACCCGTTTTGAAAATCTTGGAAACATCACGAATCGAAAGTTCCAGCATTTTCTGCATAATCTAGTTAAAACACCTGTACATGATCACACCGACTTCTCTCCCCCTTGGAGCGTAGTTATATCAACAGATCCAGGCTCAATGAAAATTAATCCATACCTCTCTTACCAGAACGCCATAAAACACTCGTCCCCATCCTAAACAACTTCCGAAAACTCTGAACCCAAACCCTTGCTAATCTTGACCTACTAGTCTTTTCAAAAGACCTCCCCCTAAAACAATTCGAAACAACACAGGGCCTCTTCCGATACAGGAATACGTAAATAAAGTAACACATTTTTTCACTAGAGATGAATTACGCGCGCCCTACACAACATAGTGACTTACTATCACCTCCGGCTACGTACTCTGGCACTCTAGATCCACGGTTCATCACAAAATCTGTCTGTAAACATAATGGCCTATCCTCATATTAAAAATTTCGCAGAAACGATGAAACTACATTCCACACACAATCTGACCGGAGCCTAACGTTAATAGAAGGCTAGACTCCACGGCCGCTTGGAATACGCCGGACAAGACCCATGCTCAGCTCTAAAGCGCTGATCACAGAGGCAGAGATATACAGTCTTACCAGCAGCGGACCTGAAATGCTATACGCTTTTAATAAATACTTTATACAGTTCCGGCGATTGGCAATCTTGCACAAACCTGCCCCCCACACTTAAACATGCGACCATGGGCTGAAAGATGATTATGGCTCATACGCTCCGAAAAAATATGAATTATCCGCTCTACAAAACTCTACGGCTGTATCCACCCACATATTCCCCCATAGGAAAGATCCTTAATGGATATCTATTCGCGTCATAAACATAAGCATCTACATGCTATAAATCTCTTTATAAGCTCATACTAGCCACACGTAACAAGACACACTCACAATTACTTCCACACACGCCCTTGATTCCAAGCAGACTCAACCATTCAAATGGTTCACTTTATCGTTTAGGTGCGGGGTATCTCGCCCGTTCCTCAATCTCTAATTCCGGACCTTATACATCCTTACAGACTCCGTATCCAGGCCCTCTACATCTCAAACTCACTAAACCACACGATCCCGAAAATAGTCCAATTCCTTCAATCCAAGCCATACCTATCTCAGGAACCTTAGTTGTCTAACCCACAATCGTACTACTTTTAGATCATACCAACTACCATTACCCTATTTGAATCACACCTAACTAATTAATTGCAACCGTGTAAATCAAAAATAAAACCCCCCAACAATAAACTAATATATGTACAATTGTAGAAGAAGGCCCACTGTCGCGCGTCTCTCCTCGCAGCCAACGGAACCCTAAAGCCTAGTCGAAGACAACTTCAACCGACATATGAAATGAAAGCAAGGCCAATAAGGGAAATATTCAGCTCCATTCTATATTACATTCATCGACATATTGATCGACAGTAAAGAACCTGACAGCTCGAACCTAATATTGCACCAAGTACCGTAAATCGAAATGCTGGAGAACAACTGCCCCTATCCAGTTACGTTGTATTCTACCAAGTGTATATCTCCGCAACATTTGGCCCCCGTTATTGATCGGCCCCCACCCTGGACACCACTTAAGAAAGCAAACTACGCGAACAGGTGGTCTTTCTTCGTTTCCGAGTCTCCGACTCATCCAGAAAATCAGCGATGTCGAACTCACGTCTCCAGTTCCTGCCATTCTCCAACTAAGCAATTGCACTCATGATGACGTTACCAATCTCATCCTTATATCCATTCAAGACGAAAGGTAATATGATTACAAAATACCGCGCGGCTGAGGTAACCCTTCGCAACTGTTACTAATCTTCGAACTTCTTAATTTGAGGGAAACCTCTCCATATGCACTCAGCCACTCACTAAGAGCATGTCCAGCAAATACTTCCAATCAAAGTTGCATCTTAACAGCCATTACAACTGGTGCCCTCTATCACAACCATGCACTGACACTATACTAGAATAAAGACTACATGGTCACCAAGACCCCTTCATGCAACCATACAACAATACGAAACCTTCGGAGAGAGGCTCCCGTAGCTTAAAACGAAACTCTTTCCATGATATCATAACTACTAGTTTCTACATTACCCGTCGACGATGGTTACGACATAGCCCAGATCCTACCACCGAATATTGTGCCAACAATTTCTATACTTCGTCTCCAATAAGACCATATGTAAACTTAAGAGCCCACGCACATCAAACTATACCAGCTAACATTACGCATGGTCTCACTCAACTTCAAAAACTACTGCGCAACTTTTAAAACACAGTCACTCCACAACATGAGCTGGGTGAGTGAATACATAGCCTAAATACCACTTTCCCCTTAACTGCTACATAGTATACCTCACTAACATACCTGTACTGTTCGTAACAAGACTTAACGAGCGGCTTTAACCATATAACGATAATCATCGTCTATATCAATACGTACAGTACCCAAACCACAATGCTGCCCATACCTCATGGCCGTCCATCTTGGTATAATTGTACACCATAGATATATCACTTCTGAAACCGCTCGTTTCGTACAGTTAACCACAACTCACAAAAATATTCACCCCATATTCCCTTATACTTAACGGTGTCCAACCCAACATGAATTCACAACTCACTTGAGCTTGAGATGCAACTGGGGGCCATGACTCACGTCGCAATTTTCGCTCAGGCAACCATCAAAACTACAAGTACAGACGCTCTAATCCACCCCCCAATATCACATCACCAGCGCAAGACTCGATTTGGTTCTAGTCGATTCAAACTTACAGCGAGCTAAGCCATTTCGATCATGAACGCGATTATAACTTTTCATTTAGCTTCCCCGTCGCTCCAACAGCCGGCTGATAACATGACCAGATAGTTCATAACAATCCATCCTCGACCGCTCACCTCAACCCACTCTACACCCGCCACAATAGCTGTCAACCGGTAAAGGACCATAATAACGCCCATTCTGTATCATCCTAAAACCTTGAGAAATTTCCTTTAGGGCTTCTTCGTATCCGTGAAGATCTCGACTCGACAGTGGCACTCCATCTAAAGGTTCCTTTTCAACACAATATGGCAAATTTACTTTTCCATCTAACGGGATACAACTCGATGTGATACTAACCAACCCAAATCTATGTTTCTAGCCCGACCTACGGCCCATACAATTGTCCATTGACCAGTCGATACCAATTTCACTACGACCAACAATAAGCTATTCATTCCCAATCATACCATCAAGACGAGATAATTCTCCCATAAGCGTCACACATCAATCATCCTCACCCAGTGGACGATCCATTCCACCCCCGCCCCAATGACTCAGACTCGCATACCGTGATGTCCAGCAACAACTTACGCCATTGCTTCTATATATACCACATTTCAAAGTAAGGATCTCGATTCCATGAGACCATGCATTGGAAACTTGAAGGACGTCCCTTCACTCATCTTCCGCAGTTGCCAGTAATCATATTAACAACATTGCGGCAAACAACAGCACCGACCATCGCACTGACGAAGCTTTACGACAAGATTATGGTGACGCTGTCTCATGTAATCCCAAACCCTAGTCGCATCCACTCAATGCTATATATAAATATCGCTACGCGTATCATATATGCGCCAGGACCAACTTACCGTGCCTTAGAGCCTATAAATACAACTCGTTTAGCCCTTATGCAGATGTCATGACGTTCGATCAACTATTATCCTCAAATCCTGTGGTATCTCCTTCATGGACAAAAGATGGAAGCTCGACAGCTAAAGCGGCCAAAAGAAATATCATCCCGGCGCTAAAAAACATTTGACTCCGTCGCAGATAACAAAGGAGTTATTGCCCTCCCGAATAGCCCCAGTCCCCTGTCAGTAGATTAAACACGAAAAACCTAGCAACTACACCTAAACAGGTGAAATCAGTCAGAACACACACCTCTGAACCAAGGAAACACTAACAGGGCAACCACGCTTAGCAACGGGTCCGACAAACATATAGATCACCTCCCTCTCTGAGCAATACGCCCAGCTCGGAAGGTGAAAAAGTATAATTACAATATAGTTACTCGCTTATAACCAAGCATAATCCTGTTCAGATAAGTAAAAATACAAACAATTCATAGTCTAGCCAGTTCGTCACCCCAGCCCACAGTAATCGACCCGGTCAGCCATCTTATCACGTTCTCCGGGCTGTAACATTCACCTTAGGACTGCCCGAAAATTCATAAGCATCAATACCACCTTGCCCCACTTTTGTTTATAAAAGATCTCTTCGATAATAGACGTCCACCTTCATCACTGCCGAAAAACATATCAGTCACATTTGCCTTATATAGCATGATGCTCTCAATCATCTCGAGTCATAGTGCAATAACCACACCTGTTCACAACACACGTCTCACCTCACAGTCATGTTCAAACAAACTCTATAATAGACTCCACACTACAGACTCTCCATACCGTGCAAAACCTACAACCAATTACCATTACGATTAACCGGCCCAGAATGTATCCTTTGAAAAACAAAGAACTACCCACTTAAAACTAAGTGAACTCATGACTATCAAATGTTCGACGATCCAACCGGCACACAGATTCCGTCTTTAGTTTCTCCGGTTTACAAGACCGAGTCTATCATTTCATTGTCCCTAATTTCCCAATTAATAAACCGTAACAACCTGAAGTTAAGAGAACCGAGCAACACAGTCAATTCTCACCGCTTTCTACCCACCTCCTTACTGCCCGACTGTCTCCCCCTTAGGGTGCACACACTAAGATACTCGGCCAGACCTAATTACCGCCTAGCCATCACACTATCCACACACACTATCCACCTCTAACCCCTAGCAAAACACTTAACGTATCCATTCCATTTAAAGTGAAACCTGATCAACCAATCGTCCAACCAGTTAATCAAAACCCTAGCCACTTACACAAAAAACACAATAAGATCCTCGAGAACCATTAAACATCCCATAGAATTTACCTAGCCATGAACTCTCTATCAAACCAATTAATGACGCATTACATTTGTTTCAAACCCACAATCCAGATCTCTTCATTACTATACACACTCGTCCCTCCACCAACTCAGAAACCAATACCTTCACTGCCACTGGAACGTTTCCTCCCTCTATACACAACAACAAATCCCGAGTTATGACTCATAATACGTGGACGATCACAGACCATCAAAACATGGCAACCACTATACGAATCCATTACAAAAATCATTTCCAACTAGATACCTAGATCAAGCTTCCAATTCTCACTTCGACACATAATAAAAACTTACTTTCCACTAAATTCGCAAGCGCAACTCCGAAAGATTCACCGCCTGAACCGTACAGAACCATGTAATGACCCAGCCCCGAGCAACCACCAAATCTAAACAGTATCCCCTTCCCACGTTTATTAGGCGGTGTTGCGCTCAAGGTTGGACACTTACCCCACCCAAAAGCATCCTGATTTATCAGTCCAGAATATACAATGTAGCAGACACCACTCAGTACGTTCTCTCACTTCACACCAGCTACCATGAGAGAACAAATCAAAGCAAATTATGATATTTCACTACAGCACACAAAAGTATTCCTAACTTCAAAATTGCAAGCGACTATTTGTCGTAAGGCTCACCTACCCCGACTAGCTCCCATATGAACTGAACGCGATCACTGCATAAACTTATTCATCTAGTATGCCAACTCCAAATATCCTGAGTAACCGTCCTCTGGCGACTGATCACTCGTAATAGTACGTACTAAAATAAATAATCTACATTGAACTTACACGAGCACAAATCCTTAATTTAACGTCTTCCAGCAAATCTAGTCGTCACCGAATCGCATTGTAACCCTAAAGAGTCCACACTATTCCACAGATTACCCCCAGTTGTGCGATGAAACGAAATTTCTAATACTATATCGAAATCCAGGTGTATACCCTCAGCTGACTCTATATATCGTCAGTCGGTCCTCCTAGACTTTACACTCCTCAACCTTCCCTAGGCCGACCCTTCAACAAAACAAGTGTTAGTAACGAATTCAATCAATAGTGTGACTCCGCACAAGTTCGCATGGTTAGCATATGCCCACAACAAATCCAAAGTTTATACATAACCAATGCTATAACATCCATCAAGACACTGTAACCCTATGATCAATCGCGGGTTTGCAACCAAAACGAAGAATGAAGTATTCTTGTGTGCCTATAATCCAAGCACCCTAAAACCTAGGTATCGACAGAGTTACAATACCAAACAGGCATGTACATTTGTTTTCGCCATAAGACCACTAGTGATAAACCTGTATACTATCCATTATTTAAGCTAAGTCGCTAACGTTAACAATCTAGAACAATTCAGACCCGGCAGGACTAAGCGGTCCCACACTATATATGCTACCCAGTACCAAACATTTCATAGACCCACTAACACCGCACGCACCGATTACCAGTTGATTAATTTACGTCACTTCATGCTCAACACATCACGCCTACATGTCACCTCACATCACCTCTGTATCTTAGGCATACCGAATCAAATCCCCCGGTCACATCTATACTACCCTCATGCCCACTATGCGTCCTCGTATTCCCTCTCACCCAAAATATCTGCCATACGTAATGATAACATTTAAGTGGACTCCGATCCGGACCATGAATTGCGGCTGCTCCAGATAAACCCGCATCTATTATAGGTACATAAATGTAAGTCCTGCAAGGCACACTAGCACGACTACCTTGATCATTATAACTACTCACAAATTCATACCTCCAGCCACTAATCTTCAAGCAGCCAACGCTCTGCGCGATCACTCGTAATATTTACACAACAAAAAAACAAGCATCCATACAAATTCAACTTGTGCTGACAGCCACTAGGTGAGATCTCGCCATGTCCGCCAATTTGTAGATTCATAATCTAACAAACTGCGATAACTGTCCTATAACGCATCACAGAAATCCTAATGCACTTCAGCTCTGCAACACAATATAGAATAAACAACCCTACGAAACTCGGCCTGGACTATCACCGCGTTCCATGATTCGAAACGTCCCCATCTCGAAGCATCCTCTAAGCGCCTACCACCCCGTCAGATAGCCCTCCAGCGGGCCCCATTCCTGACCCGCAGCTCACACAAACTTAAAATACCCACAAACTGACCAGAATCCTTCGTCAACCTGAAAGAACTCTTATGAATGTTTGCGCACTAACCATAAAAAAATGAGTACACAGTGTTATAAGTTGCATTTCATTCAGCAATATTTATTCTACTTAGTCATTTCATTACAACCGTGTGCCTGTGGAGCAAATATCTAAATTTATTAAGGAAACACCCTCATATAGGATCATGCAACTCACTTAGCATTACAACACCACTACGATTTTGATTTCCTAGAATTACTCACAAACTACAACACAATAGAAGCCACAGTATTTTAGCGCGTGGTTGATAATTACAGTAACCACACCCGATAATTTACCTAGGTATATCGCCCAAGACCGTCCGGCACTACAATTATAGGTCCGACTACGCTAGATTGTATCGACGCCAGGTAAAATTGCGATGTACATATCACGCCCCAACCGTTTTCATCACCCAGCAATAAGACCGAACTCACAACTTAACAACAAAATCTAATAGAAATAATCAATACATCTTACCATCCATACTCTGATACATCTAATCTCATAGGAGTCCATACACAAGCCTAACAATACTCCCCGTGGCTTGTTCACACTGAGTATGCACATAACCTTTTCGTGGTCTTCTGAACCGCACCGAAAGCAACTCCAGATTAGGCTAGAAATCTCCTTTCATATATATTATAACAAGGAACAATGCACGCCGCACCTATCAGCGTTTAGCCACCCCGAAGGGCTAAGAACTCACCTGATATCCGTGAAACACTCAACCTACTTGTACCCCTTAACAAACACGAAATTTCTCCTACCGCTTCAATTATTCACTAATTCACAATCCTTCCACTATCATTGTACGCCAAGCCCCCAGAATAAAAGATAAAGTTTCTTTCACCTCACCTAGACAGTATATGAACGTAAACATCCGATCTTCGCGCTAAACAACCTCACCGACCAATAAAGCGGTTCACTATGTATACTCCCGACACTTAACCTAACTTGATAGCCTTTGTCCGCAGAAAAACTCTGCCTAGTCAAATACTAATATAAATCAACCCAATCAAAGGGCCCATACTACCACGCACCGACTACGACTCAACAACAAACTACCCTAAATCACCAAGATACCCCTTTATGCTGTCTCCCCAGAGCGCCACCACCATAAACGCACACCCTAACGGCTCCTTAGCAGCCCCATATGTGTCTATTATATTTATACAGGGGTACCTAATGAGAAACCATCTACATTCGAACTTGTCCCCAACGCCTGCGCATATCGGGTAGTCCATAACTCCCCCCCTAAAAAAACTATTAAATCGCATTAGAGGCAGCAAAAGAGTTACTGACCCCCTAAGCCCCGGGGCATACGAATTAGGCCCGTCTCAACCCATTCATTGGAACTACATAGACAGACACCTAATTATGTAACGAAAAGCGCCCCCCACACAGACGTCAAACTATGAAACAGGCAATCACACAAAAGCTACACTCTTGCTATTACATAGCATGCGGATAACTAACACCGATAGTACGAACGAAACAACGCCATGCACGTAACCCATTCTCAGTCCACCCCCATACTTCTAAAACTAAACAAGTATTTGAATAAACCTAAGACAAGTCGGAATCCTCTTGCATAACATTCTTAATCAAGGCTACAGTCATTGCACACCACGAAAAAAGACTTACTTCCTATTATATGTCCTTAACACACTTAACTCGTGCCCCTTTTAGACACCCCCCTACGCAAATTCTAAAAATTATCCATGATTCGACCTGCAAAGATAACATCACCATAAACTAATCAATCTCTTATCTAAAACCTGTCCCCCACCACCCTCACCAGCATGTACACGCATCGCAATACTATTAGGAATCTCCCCCGTCATTTCTCTACATACTTACAAGTGCTAATATTAAACACATGCAATCCCCTCTTAAGATCAGCCTTAGAAGAGTGAATCCCTTAAGACCCAACCCCCTTGCTGGAATCAACTACGGAGTGGAATTTAACTACTCAACTCCTATATGGAATTCCAGATCTGGCAATACCGCGTCCTGCCTCTGTGTCCTTCTATAACAGCACCCCTGCCATCACGCCTTATGTCTAGCCGAGTACTGACTTCCGCTCCTGTCAATACCTTATCACAAACAAAATTGACCATACCGTTGAGGCTAGCCGCAAACATACACCGATCCTGCCCAATTCTAACGTTTACCGTTGACGTGTTCGTTACATGAACAGCCAGAACACGTCATACTAGTACATATAGATCCCGACACTCGGCCAAACTTTTACCACCTGCGTAAACCTATGAGATCGTCAAAATACAAAATAAATAACCACCCGATATACACTAGTACACCAAGCCTCCGCTTACCTAGGAACCACACCTGTACATCACAAATGATAGGCAAATACATACGCCCGCACTATTCACACAACAAAACCAAGTAGGAGAAAAATCCCTTCCAGATCTTTACAACTTAGAAATACTCTCATTTGCTTGTTAATTGATTTATGATTAAGACAGTTAGCCTGTTAACAGAACCTTCCACTTATTAAGTGCTTACGGGGGTCACACACATAAATTTAAACTTTATACGTGTTAAAAAGTGTCCAATTAGATTCGCCAACACTTGCAAAGACCAATTTGCTAAAATAGCCTCAACAAGAAATTCTAACCTTTAGCATTATACTATCGAATACTTCACCTGATCCTTCATAAGCGAACATGCACTCCCCCTTTGTGGCAGGCGTCTCTCGGATACCATTTGATACAACACACATCAACCCGCCACTTAAACACACATTAACATAAATCCTGCAGGTATTCTATGTACTATCGCCGACAAAACGTATGTGTAAACACATATTCCATCATAAGAATCAAAAAGAACTATCGTTAGTATCAAAATTTCGTAGCACATCGATATTTGGCCCCTTATCCCCATCGTTGTGCCGGGATCTACAAAGCTTCACTCAACCACTCCCTCCTTACTGCACACCCGTACGTTTTCGACACTTAATGATGATATCGCACCACCCCTAGCCGTCGCCCGTTTCCATTTGTGGCATTTTGTAACTGCTTCGAGACCTGCAAAATCTTAGGACGGTAAAGACATACACGAAGGTGCCAGAGACCACACGAAAACAACTAGGTAAATGTATGTCCTCCAACTGAAGTCAAAAATATGTCATCATTTGGTTCCAACCGTATTACCCAAGGGTTACACCCTAGTACTCTTACGATAGCACAGCGTATTCAACCTCGCCTAATAGTGACTTCAAACTTTAAAAGACGGTTACAGGCGCAGGATCTATTAAAAGCTAACTGACTGGTCAGCAAAGAATCATAACCTATCAGTAACTGTAAAAGGGCACGCTCCTCTACGATAGGTCTGTTAGCAATCTAAATTATACAACAAGCTGCATTCACGAGAGCCTATTAACTTAAGCATCCAATGTTATACCTCGTGACAATGACTTGCTAACGCCCCTAATATCAGAGTGAAACTATTTAAACGCAAACTCATTACATTCATGACTATATGCAAGTGAAGTGCTACTTACCTTCTATAGAACACTACTCATATAAAAGCAAATATTCACGGAATTTCAATTTCATCCCACTCGCGTCAGTTACGCTTCCATCTTCTCTAAGCATACCACCAACCGGTCCCTCAGTACCCCGACGTGTGCCTTACTTCCACTCAACACCCCTTGTCCTCTATTGCCCTAACTACGCTCCACCACAAACGGCACATATGCTACAACCCGTGACACAAAAAATCTACCCGACCTCCACGACTAACAAATCCAAGACCTCTTCACAATAGGTCTCTTATCGAGCCAAAATTTTCCACTGTTCAATCTTCAACAGAAGGAGAACTTTCAAGTTGAATACCTAAAAGTTATCGACCAATGACCTATTTGCTTTGATCTCATAAGACAAGCTAATCTGCCCAAGTTGCATTCATTATTGCAGCTAATCACTAACAACCTCTTTACTCCCTTACCGGCACCTCGGAAACACTATCCGCTCACGTCGACCGACATCAACCCAAGTCTATAAATACTCATTATATACTCTAAACGTGAATAAAACTTCAGTTGAGGCAAAACTGTTTTTAAATAGCCTACTACTAAATATCGCGCACTTCTCATTTACACCCGCACTTCATCCATGGTTGACCTATACTACATCGCTGCTCGCATGCCTCTTTCTATACCTTCTCGCTCTTACCGGCCTGCTCCCAATTATAAGATTCAGACTGCGATACAATTTGTAAAACTGAGATAGTCAATGCGACGCTGAGAGACATTACACGTTTAGACAAAAGGGAAACACAGCTTTTATAAACAACTACGAGTCTTCGCTATTCTGTCACAGCCCTTTATAACAAATGACTAAGTCATAAGGCTCCATATTTATACTCGAGACTGCCCTCCGACCTGCCGTTGTGCGAAGTAGCTATGCACTAACAGAAGCGGAACTCCCAAGTTTGTCACACCTATACTATATTTCCCACGACCACGGATACACTTCTGTCACCAATTAATCCATCCAGATACTGCTCCATATTCCCCGATTTAGACTAACACATTGCTCAAACTTGAATAGGAAAATTATAATCGAGTGAATCTAACAGGTAACTTATTTATTTCCGCAGGACCCCGCATCTCAACTACTCCCGACTACCATTGATCAAAACCGCGTCCACTACATTGATATGAATCTGCAGAAGAACCCCTAATCACAAATAATCGGTAAAACATAACCCATTAAGTGACACCATAATCCGTTACCCCACATAAGCACCCACAGCAGAGATTATTGGCCCGGTGCGCCTCCCCCTACTCCATAAGTCTTTAATCTGTCCCTAATCGACCAAACTTCCCGTTCTGAAGGCACATCCCGGGATCGAAATCACAACTCAACACACGTTAGTTTGACCCGACTAAAGATGGTACAAGCTGTACTCTCCTTCATTCAAACTTTCCCAAAACACCGAACCTAAATCAAAGTTTCTATGAATTTTCTCACCCCAACTAACAATCATCGCAAAATAAATCTCCAAAAAATTCCCTATTCTGATGCTTGATGCACTAGTATTGATAAGAATCGATCTAACACCAACACGTAATACTTTTATCTGACAAAGAAATACGTACCTCCTAAAACCCGCACAATTCTATAGCGTTCGGTCACAGTAAATTTGCATTGCCTCTCGTCACTGCCCAACAGACGGCCGCTACCTCAGCGGTTAAGCCACCCTCCCTACTACGAAATCTTCACGCAGACCGCCCCGAACCCTAAGAACACGAAGTCTAGCGGTATCGATGTAGCGCTCTCGAAAGACGAGCCCCGGATTGTCCCAGCAAGTGACAAACAAGTAAAATTTATAAAGCAAAACCCTAATTACGTTCCCCAACCAACGCTACACTTTTTAGGTTCAGGACCCTACAGGGCACACACACTTGGCATCGATAATTGATGAAAGTAACTCGTAAAATCCACATTCCCTTACCCCACGACCAAAAATAACCAGCGAATCTCATACCCGTCCCTTCTATCGTCTAAAACTGTAAAAGATCGCATTCATCCCACAGAGTTGGTCGCGATCGTGCCCTATAGCAGCGTATTAATTGACCAGATCGGACCCCCCAGCAAGCTATTGCGCCATAGTAGACACAAACTACAGCAGTTAGTAACTGTGAAATGGACTTAACTGCAACGCCAAACCATATGGTCACTCAACCATCACCAAACGCCTTGCGATAACTACATTTCAGCAACACTTAACAGCATGCACTACAATCCAAATGAGCGCCCTATGTCATAACAAATTTACCTTTTCTAACCGTCCAGGAGAGTACCATATCAGAACTAATCCTCTCCGGCACTAATTACAGTCTAGGCACACCATTTATCCTTGATAATAGCCCCCTATAGGTCACGGTCTAAAAGTGCCCAACACTCTTACAATTGCAACGAGCTCAGATCGTCATCGCAAATTCACCAAAAAACAACAAAACATTTCTGCCAGCGAACGGATCTAGGTGTGCTCTTAATTTTCCGAGATCAGCAAGTGACCGCAACGTGATCCCATAACAGACTGCACCCCCTATGGTCAGCCATTTCCGCTAGAAAATCAAAATATCCCACGAAACTCTATACATTTAGAGCAATCCTCTCTACCGGTCTCTAGTTAAACAAAGGGGTCAAATAAACACAAATGCATTCTTCCAAAAATATCCGTTCTAGCGAGTAAATAACCCACGCTTACGCTGAACACTACAGACTGCATCATTAAGCTACCATCATTAACAATCAGCGAATATACTCAAAACCCAAAGCACTTGGACCAATTATACAGACCCCTAAATTACGATCACTGTCATAAATCACTAGTCCATGGTACCCTATCCGCATTAAGCACCCGTCAGATGCACCACTCAAGACTCTACTGATTCAGACGCACCATTTTCAGTGCGACATAAAAGCAACAACCCCACATCTTACGAATGAATAAAACAGCACAATCATACATCATTTTATGCCTACCTCACGGCTCTACCACGAATCGCTACATCCCCACGCTCCCAAGATGTAATAACGGGGCAAGAACAAACAGGAGATTTGAATACCCACACATGGTCTGATGTAATTCAGGGAGTTGCGCTAAACCGCAATCCGCAACCAATAAAAAGGCCTCAAGTTCCCACAGTATGACTCAAGGAATACCTGCAGTGCACCTATGAGAATCTATCTAAGGTACATCTCCTGTTGCTACCTAAAAACAATTACAGTCACACCACACCCCGACAACGATGATTTCTAAATCCAAGACGCCTTAACAGAACTCTGAGCTTAACCCTTAATATGGGCATGAATAATCGCCTTTAAACAGAAACAAATTCGCAAATGATTAAATAACAAGATAATTGCCTATCTTTGGACTAAAGGCCATCGTGTGGGCACAGCATAAGCATCAGATAAAATGTTACTAATTCCAATAGCGGCTTATACCTCTGCCTCTAATCCGCCTCTGCCCAGATAAATTACCAAGAGGAATTCCCATTTACAATGGAAAGGCCCACAACAGTTGTAAAACTCAAGACCCGTAAGGACGCTGATCGCATGTCAAACCCCGCACCCTAATTTGTTTAGCCTTAACACACTACCCATCTTAATACAGTAAGGATGCAGATAGAGATCTCATTGAGCGGTCAAATGGTCACTGAGACAATCATTCTAGCTATAGATTACAAATAATACTCTCCGGTAATTTTGTGGTTTTCTTAAGCACCCTGAAGATGCAAACAACATGTGGAGGTCCGCCGGTATACACCCGAGCCAATTCCAATCATTCCATTGTGCTATCCCCAATAAGCCACTTCGAGACATAACGTGAATGTACATAAAGACAAGGACTCCCTCATCAGTTTTACAGCAGTTCCCAACCTTATCACTATATCCACTACTCTAACCGCGAGCATATCCTCCCTTATAAGTCCCGTACCGAGGCGCTATCGTCTTAAAAAAAACGCAAACACCCAACGCGCCTACCTTATTCATTATCAATCGGATCAAATGTTTTACATAGCGATGGCCCTACTGAAGAACTGCCTAACAATGCCCCAAAATACCACATATTCTTTTTCAGCGGCTACGTCATATCCACAAAGAAACTCTTGTACACCAAAATTCTAACTCCTAATTTCTTTCCAGCGATTAGCTATCTGAACTGAGCCATACAAAGAATTCGCCTAAATCATGTATCACACGACGCCTTATAATCAAATCATCATGGGCTATCTAACCGGTCCCATACCGTGTCCCCACTACCGGCACATCTATAAGCAAGCCAGCACAATTCTAACCTCTTGGCTGAATATCAGAACCCAATTTTACCAGTCGCTATAAATTATTTCGCAACAATTCGACATACATCCTAAAGCTGTCCTCTCTAACAAATCCGTCGAACTGTCGCTAGCGAATACCATTTTTTCCGTGAACCGGAGTACTACCACAAACATAGACCAACCTTTATCATCAACAGACTACAAAAATTTAAAATCCACTAGAACTAAGCCGTCCTACTCCATTGCCAATAGACATACCCATACGCAACTCCTGTCCTTATGGAACACTCACCACCAGTGAACCATTATTAGAAATTACCCGAAACGGCACTGATCAATAACACAGCTACAATCTCACGGATGATACGACACTGTTTCCTCATACAACGAAGATTCAATACAAACCTTTATTTTCGACCTCCCTCATATACGATACCAGTCGAGTATCGGAACAACAATAACTGTATACTAAACATATTGCCGATGATACGTCCAATTACAAACACATCCCCCTCCAGAAAAATATTACTTAACCAAACACACGAACTACGAGTCACAGCTCCCGCCACCCTTATTTTCAAAACAAATAACACTAAGAGAAGTGTGCTCTTATACATATGCGACGACACGACTATTCACCGTATAAATTACATTGCCTTTACTGTCGTCACTACATCATCAATTTAAATTAGCCAAATCAACAAACACAATAAATTACTCAACAAGATCCAGCACGTCCTCCAGTCGCCTACCATCAACATTCGCCATAAAATGATGTCTTGCGTAACCTCGCGTTTCACTCACC

hippo_random = hippo(randperm(length(hippo)))

hippo_random =

AGCACTAGTGGGAGGTTATCTAGTTATAAAAACATATCTAATAATCCCCGACAACCGACGGTCATTATCTACACCCTCCAGTAAAAAATGACCGACATCACCCTATTCCGACTTTCCTATCGGCTAGACCATCCTGTATATGAGTCCACTTAGATCAACGTTTCATCTAATTTTTAAATAAAAAAACCCCAATCATCACCCGATAAAACCCAATTCTTCCAATTTAACCTACCGCAACCCGCCTATACCTAGAAATTATACCGTACTTAAGAAAATAACACAAAATAGTTCTAAAGAAATCGACTTAATATTACCACAGCTCATATCCCTGATGTCTTATCGTTACACCCCGTCTTAAGTGCGACAAAACCTCGCTCTATATACGTATCCTCTGTAGTCAAACGAAAGAATCGTGTAACCCCCCAAAATAGACCTGCGCCAGGAGACAGTGTTGGATCACAATGAGCTTGTCCAAATTGGTTACTCGAACTCTTTTCCAGGCCAGCAACCTTCTATCATTTTCATAGCTGTCATCACAAACAATCGCCTATTAAACATAAACATTCGTCAAACATAATTACAACACGACTTGTACCCCTTAACAATTACCCAAGACGCACGACATACAATCCACGGAAATAGTTACCGCTAAACCAATGCAATATGTAGTGCAAATCACTGCTTTTATCTTTCAAACACTATACAGCAATAAATAAATTATACAAGCGATATTACGGCGGCTCACGCCAGCACTCAATCCGTGTTGACAAATGCACCATAGTCATACCTTATGTACAATACTCGCCATACGCCCCAAGCGAATTTTAGGTGTCAAGGTCAACTAACACCGCCGTGTTCCATTCCTGTAGGTATGTTCGCTGCAGCACACAACGCAAATGGTCAGCTACTCAGTGTACATCGTCCATTCACACCACCCAGCTAAAGCCTCTACTACTAGCCCTCCGCCTCCAAGATCTACCATCAGATGCAACTAACTATATAGCAGATCTATATTATTTTCAACAGAAATCCCAAGAACCTACCAGTTATTACAACCCCACTAATTGGATGGCATAACTTCCGCTAAGAGCATGACACGCTCCCAATTTGAGTACCAATGTGTCCAATACCGTCGCACCTACTCAACTAAATAAATGCACACATTCATAGGCACCGATAAAAAAACTAGAGATCACATGACTTCTTTTACGATATCCTGAACTTATCGGTCAATGCGCTCATATCCCACTCATTGTCATCGTCATACACCATCTCTTCTACGAGGAAATTTACCTACCCTTTCCATTGAACAAGATGTCACAATATATGGTATTCCCTGAAACCAACTGGTCACTATCGATTGTGGAAAAACTATTCTAAATTTCTCGCATTGTAATCCACTTAGCCCCCAGTATCCCGGAATCAGCCTCTCGTACTGGGATCATCACCCAAATTTTTAACGACCCTTCTAATTAATACTAATCTTTACAACAGAAGCTCTGCCAATAGCTTAGATCAGAACAGCAAGCTCTTAGTGCTTGCGCAATCTCATCTTGCGCAAATTGGTAAAAAATTAAATACCTAACTTTATCGAACTACAGGAAGTATCAATGTCCACTAGAAATCTGCATATATCTAACAAGGTAAAATGCAACTTTATGTAACCCTTTTCTTTCGGTGTCCCACTCCTCCTGCTAACATATGAGCCTGATACTAGACGAACAGACATGTACACTGACATCAGGTAGAATCGGCTGACGCCATAGATTAGACCTGGGTTAGTTCTGGCCATAACATAGACCAGTCGAACGCCTATTACAATATACCCCTTCTACAACCACGATTGCCTATCTACCATCACCTTCGGAGAAAATATGCTTACCAAGGCCATACTAACCCAAGGTTGATTACTCAAACTTCCCACTAGTTCCAAGCATCGCCCATGCCTAGATTATTTCCCTGTCACCGCGCACACCTATCCGTACATGATACAGTTCCAGGCACACGCATAGGTCCATGAGAAGTTTAGTCCCACATTGGCTATATTACTTTGATAACTATAAAACCCCTGACGTGAACGTTGACCACAATTACCCACTCGACAAACTGTACAACCACGTTTTTCTTAATCCGCACTTTCGTTATCCATGAGCTCTTCAAGTAAATTACCACCTCAACGCTGATGAACATCCATCGCTTTCCAGTCACGTCTCCATCCTTCAAACCGCACGCAAGGCTACAACGCCGTCACAACCAACCCTATCAATTTAAATCATCTTGGCGTCTTACTCCTATCCTTATCACTCCTTAAAATAGTACCAATGATCACGTCACGCGAGTATCTAGTACACTAATCATGTTCAAGAATGCTCTTTCTGTAAGGCAGTGAATTATTATATTTACAGGCAGTTAATATCGTTTGCGTATACCCTTATACTGTACCACCACAGATTCCATCTCCTACACCTTTCACAATCACGGTCCAGCTACAATATGAATACTCGTTGAAACATCCGGTACTCATCTTGACCACGACATACGGCACGTTGGTTGGAATAAAGGTTTTTTCGTCCGCTAAGTTACCCAAACAACGTCGGACTAACCTTCTTAATTTGCGCCCGACCCGGAAGTTGATTGTCTCCTTATATAGCAAAACCCCTATTCACTCGACATGAAACCTTAATAAATACCCTATAAGTCCCAGATTCCGGCCGCAATCTGTCAGGGCCGTCTCTGAAACAGCTCTCCACGTTAAATACGGGCTCAACTTACCTCGCCAATATCATTTTCTACCGACAGGCAAGGCAAGTCCCATAAGCCCGCCGCTCGTCTCTAGTCTTCGACCTGATGGATCTCAGAAATTTTTTCTTGGATCACTCAGATGCTCACCTAATGCACAGAATATAGCAGTACAATCAGATTAAATAATCCTTTCCTACCACCATATACATCAGATTACCAAGGTATCCTCCCTCACATTCAAAGTCTTAACCAAACCACAGAACTAAGCACCCGATACACTTTTTACCACGCCAGAACAACTTAGCGAGTTGCCTGACGACAGACCAAACGGCTACCTTAAGACAAAGGACCACTACGCACCACCCAAATTCCAGTGCAAAATGCCAAAAATTACAAACGATTTAAAAAACAGTAAGACGTCAAAACTGTCCTAAACCGTGTGTGACCAAGCCAATAATGAATCACAATTACGCGCGTATTCCAGTTCGAATTCTTACTTCGACTCCCGCTCCCCTGGCTCGACAATGCCGATCCGTGTCATCTCCTAAAGTATTCTAAAACCACACTGTAGACAAGGCCTATACACCACCCCGCCAGAACACGACAAGCACGGTCATAACCAAAAAAGAAGCAGAAGAATACAAAAGACCATCGGCCGTACATGAACCCAGACCTATAACACTCATCCATTAACCACAACATAAGATAGACCGCCATACTAAACAATTTCGACGTTTATGGCGACCTAAACTCCCGACTGTCCCATTCGACCACACACTAAGTACATCCTCCTCCCGTCAACACCGGCGGTACGTACGGAAAAACGAGCCGATTCAGTTGCCAACAAAGCGAGAGCGCTGGTCGGAATCCTCCGGTCAACGTACATGGAGTATCCACCCCAATTACAGGATCCACACTCCCTGATTTTCCGGTACAAACCACCATAACAGTCATCCGAAGATTTACGTCCATTACGTACAGTCTAATCTAACGAAACGATCACATATAATCCATATCCAACCAGAGGGTCAAACGCATGATCAACCTCACTTCCCAGTCTGACCAAAAACATGAGCTAAACACCTCTTGAAAACCCAAATATCGTTCATTCGAGCTGCTATATCAGGCAGATTTCGCCTAAAGGATGTAACATTCCCGAAGCGACATTTTCAAGTGATCAAAAAGACTATCATTAACAATCGTATACGTAAGACAATACGATAATCCTATCAGTGCCCCTCTCTTCGTCTAACTATCACAAGCGCCAGCAGTCGCATGAACCGCGTGTTCCCTATGTCAAACTAAAAATAACCCTTAGCTAAAATTAACATCAAGTCTAGTTCCACCCTTAACATAAATGTCGATCAACTATTCAAAGCTCTAGTGATTTTTCCCTGTAGTAAAGCTCATACCTAACCCAACCATACTGTCATCTGCAAAAGTCGAACGTCTGGGGCAATCAATCAGAGACTACCATCGACTCGAGCAAGTACACAATACATGGGAATAAACGACACAATCAATGATACATAACCCTCATGACTAAGAATCTTTATCCCCATCACAAACGGATCACTTCCAAACTAAGAAACCTGGATCAAAACTCACTTCATTACGCAACGGACCATTCTCAAACAAGCTTATCAAAGACGACAGCAATTAAACAATGTCATAAAAGCCTCCATAGACTATATCATATAAAAATCCCATCCCAACTGAGCCGATTCCATGAAACTGCTTGAGTGAATAAAAATGAAACCAACATATGCCACCATTCCAACGCAGGCGACCAACAAATTCAATAATAAAAGCCTTCATCCAAGACTACTCTCCTATCTATCCCACATTACCGAATACTCCCATGTTATAAACTGTAACATAACAGTGAAGTTTACAAACAACGAGCAATAGTTCCGAAACTAAAACAAATTCATCTACCAGGACCATAACACTCAGACCTGTTCCCGTGAGTGGATCCTAACCTGATGCTTTTAAAGAACCAATTTACCCTCCACATCCAACCCGGACTGACGTACGCATTGAATGGCTCGCTCGTAAAACCACTGTCCCCCAGATACGGTAGAGATAGTCAAGTTATCCACTTAAACATATGACCGAATTTACCATTAATGCTCCCCACCAGCTTCTTCAGCAACGTTGATTCACTTACAAAGACAAGCACAAAAAAGCTCTAAAACCCGGCACTTACGCCAGACCTTTGTTGAAAACGTACCATCAACATCTAGACGTGAACAAGAGTATGACCGATACCCGTAACTACACCCGTAACAACTACAGCAAGCAACTTACTACCTACCCCCATTGAAGCTCAAGTGTTATTGATGTTGCACGCCACTCATCCATTAAATCTATAGCCCAGACAACCTTTGCAGAATTATCCAAGATGACATGATTCTAATAGTACTGCCCCGCCCGCCGGGGGAGAAACAGATCCATCAATACACAAGAACAGACCTGAATCAAACACCAAAAAAAGCATTTCCGATTTGGAAAAGCAACAACCTTCATAAACTAATGCAGAAAAACCCCCAAACATCACAGCATCACAGGAACTTACTAGTGTCAAAACTATTCTCCCTATGATCATGGAAAAAACTAAGTGCCACTTTTGCTTCTAGCTAAAAAAAGAATACTTGAAAAACGCCCCTTGGTTGTTAATAAATGACTCTATAGCCCTAAAGGGTTTGTCTACACTACCATCTTTCACCACATTCCATTCATGTAAAACTATCAGCTCACATATCTGGTAACGCTAACTTTAAAGTGTTGGCTACTCTCTACACGAACTCTAATAACCCCCCTTCTAAAACTCAATCCTACCACCAACGGCATGGCATACATAAATAAGATGAGTTATCACTATGTAATAATATACCTAGCACCAATCTTCACGCAGCCCACCAAATACTCCCAAAACGTTAAGATATTAAATCACCGCTTTTGAGATTCTCTTCAGTTGGAACAACATTACCAGTTTCCCGCACAGCATACCTCATACCTACGATCTAGCAAAGCTATCGCCAGGAGAAAATCATACACGTCCATGAACATGCCATTGAGTACTTGGAGAAGATAACCCCAAACCTCTTACCGTCTTCCCTTAGTAAAACGACCCCTGACGAGATCCTACAGACTGTACATTTTATTAAAAGGTGATCCGACAACAACTCACCAAAAACCACCCCTGCCCCCCCATGATCCCCATTTTGGCCATCGTTTTTCTACTCACCTACGAAAACTACCCCAGTTAACAGCCCAAGTCTTTGGCCGCTTTATTTGTTAAGTGCACAAATCGTCACTCTGAACAGTCACCCATAATATCGAATAAATACATCCCTCTACACAGAGGAAACAGTCCATCTAGCTGAGATGACGCATATTTTGACCTGTACCCGTTAAATAAAACAACAAAGATGTGTTCGGAATATATTCCTTTTAAACCACGATACCGAACCCAACGGAATAGTAGATTTCGACATATCTAACATCCAAATTAAAAATAGGCAATGTCACTCTTCCTAACCAGTGAAGACGTAACTTAACCTTGATTTCGACATTAACTCATCATTACTCCCCCAACATCGACATACCAACACGCTTACATTTCAGTTCAACGAATATACACAGCAACTCAGACTCAGCAGTTGACATACCACCCTGACTAGACACCAGTTGCGACCCAGCTAAAGCGCATATTGACCCACTTTCCCCCAAAATAGCCAGCATTCTGCGTCAAACTCATATAACCAGACCCAATTCACGCGCCTCTGAAGGCACAATAATTCCACCATTGGACAACCCGGATATGAGAACTCTTATCCCTTCCACTACCCCCTTAAACACTTTAGAGATTCAGCCAACCGAGATTTCACAGGAAATCAGCCCGAGAAATCTGCCCCTAATAAGCTATGATGCCGAATGCTTAACCTTCCCACATCAAACAATCAGAATCGCACACAACCCTCCTCCATCTCCCTAACGCTCGTCAAACTCAAACTGTCATGGCATAGCAATCTAACACCAGGCACAATCTACCCAGTATGTCGCTATCTTTTTGTCAGTAACTAAAACAAAAGATGACAACCTTAACTTTAACTTAACCCGAAGACTTACTCTTCTCATAAATCGGTTAAAAGAACCTAGTTCCGATATACCGGCCACTCATCAGACGGTCTTGTTTAGCTCCGCCAACTCAATACCATTACTCTTCACACTCTAGCGATCTACGTTACGATCGCGAATGTTCATCCTAAGTTGATCACATCACAAGGCTTTACACTCTATTACTCACGGTCATCTAAAATCCATCTGCTACTATTCTCTGTCACACTCGAAGCTAAATCTCCTGTTGTACTTTCAACCTCTATGATTGGGGAACCGCAGACACAACACCATTATATACCCAATCCTTGCATCCTCCAGCAACCAAATCAACCGTAAATGTCCACCCAGCTTACCTGAAAAGCCCACCGTTCAGCTTAAGTAAACCCTTGACACTTATAATTACCCACCAACCCGCCAGTCACTAATTTGTTAATTTTTCGTGTACATGTCATATCAAGGACACCATTTGCCCGTTCTTAAGATAATAGCAAAAATTGCAAACCGTTCATCAGTACATCTGAATCAATCCGTCCCTCCTCTGCACATGACAATAAGATATGGTAGCTATACTGCCATAGTCACTTAGGATACATAAACTTATCAAATGCCACATTTTCCATCACCTCACTCCGCATAATCAAGTTAAACCTAGAAAACTAAACACGAGAAGTACACACATTGGAGACAATCCCCCCCAAATCAAATTCGTAACAAACCACCGTATTCTCAACACCTTGTACATCGTTTTTGTCATCAAACGACTCCAGTGAGTAATTAGATCTAACCGATTCTCTGCAGTCCTCAACCGGTGACATCGGATTGAAACCGTAAACCGCCTTCCTATATCCCCATATACGGTTACCAGATAGCGATTCACAACTTATAGCCGCCACGCCAAATATATACAACCCATCTCTCATTAATCTAATCATCGCGTTCCCTACCTCTTCCACAAATGTCAGTAAATTACTACAGACGGGCTGGCCCAAAATCCAGCACGCCCTAGGCACCCGATCGTTCATGAAACGAATAGATAAGGTGTTCAAACGTCCTAACAAGTCCTCCAAAGAATGACGCATCAAAATAATTGGCTCCGCAATTCACAACATCAACCCACCGCCCACACCCAGCAGTCTTGCAAATTTACCCGAAACATAAAAACGATAGGTCCCTTCCTGTCACCAGGAACTCCGCAAAATTTAAATGAATCACTAGTATGATGCCATTATATGCTGTCACGCCCCTCCCTGCAAACCCATACTGCCTAAGTCAACACTACCAGCTTAGGTCAAACGCACTATCCGCCAAAACCGTCCGTCCTGCATGTGTCCTTACCTCGCACACAAAATGAATCCAATTCAGACTCTGGCTAACAGCAATCCTGACACACAGCAATGACCAACTGAAATGAGCCCGGACATATCAATACGTACCCTTCGAAACTAGAACTCAACCGTATAGACTCTACTACCTTTCGCATATAGTTGGACTCAAGTCTGATACTGCCCACTCAGCATTTCCGAGTCTAACCGCAGCAAAATACAATCGCCCATCAGGCTAGCTGATCCACACATCAATTTTAATACACGTTGGTAACTTATAATTGATACTTCTAAAATCCCTCCGTATACCACCCATTACAAGACCACCAATTAGTCCCCATATGTAGCCATACAATACGACCTTATAGACCTCCCTATTTTACCAAAAATACATTACTCTACAAAACCTACTTAAGCTCGCCTAACATAATTCCCTTGGGCCAGAACGTTCATCCACCCTATATCACCACAGTTCCGACAACAGCTATCTAACCAACAACAGCAATAAGTTGCTCTTCGGCGATCCAAACTCCAGGTTACAAAGTTACCGAATGTCACTGTAATGCGCTTACGTTGAGAAAACAGTTGTACTGATATAGTGCCAAACTCCAAAGTGCAAACTTAGGGCCGTAGCGAGCCCCCTTAAAAGATATATAAATCGTTATATAAGGAACTAACGCTCGGCGCACAACGACATCCGCCGATGGGAATCTGCTCCCCGAAGAACTCCGATTCTATTTCAAACGCATTGACCTTATACCCCTAGTTTATACCGATTGATCAAAGCTAATCGTATTTGAACCAACTCTAATTACACTCATAATGCTTCAGGCCGCAGGACAGCTCCGCGAATCCACATACTATATATGTATTGAGCCAAAACCAATCATAAAATTCACCTAGCCTCAGTAAGTATTCACAGTACTAAAACTAATTATAGCTTGGCAACAGCTCCTATAATGAGGTTTGACGTGGACCTAACCTCATTCAGGTGTTCCCTACAAAAAAGACATTAAATCACCTACACCATCATTCTCCACTAACCTAGACTTTTATCATACTACTGACACCTGAACACGTTGCTTACCCACAAAACATCAGAACCTTTCCTCACCTAAGCATCAACGCAAATACTAACTGGGATCCCTCTGTGTTGTTGACAAAACGTTATTTACAATGCTCAACTCATCATATCCAAAAATGCCATCATATCAATCCAACGGACACTGATCGTCGACAGATAGTCATCATCTCCCCATGCCCTTGGTATTACGAACAACGCAAGGCAACATTAGCACGGAAATTAGAACTACATGTGGACAAGGTTCCTTTAAAATTAACCTCAAACAATTATAAACCTCAGCTTGAACCTCCTACCCCCAAACAACACCTGCTACCTTTGCCACCTGATGCGCAAGTCCCATCCAACATCTTCCTTTATTATCCAGATACGTTCTAACAAAAATGACAATCCTCAGGTCTAAAACTTTACATTTCTAAACACTCTACCATACCCACCGCCTAATAAGCGTACAAAAACAGCGATAGTTATGACATCGAGTGCTCGCTACTTAGTCTGTGCCTAGTTTCGTATCAGTCCAGGTTGCCTGTTACACATTTCGCAACACGAGCGAAGCTCAGACTCAGTACATTCGTCATATTTACACATTACACTCTGGTTTCACTTATCAGCTCACCGAGGATAAAAACTTAAACACAAACACGGCAAGTTTCGCCGATCCCAGAGCCTGTAAAACAGGAGCGAAGAAACTGATATTAAATCTGCCTACCACCCAACCGTTGATACAATCCAAATTTCAGACTGTGTGTCTTTGCCTAATCATCTACTTCAACGTAATTTGCTACCCGTCCATAGACGTCTAATACACGATCCTTATCGAAACCGAAATTAACCTAAATGGTAAACCACTAACTCTACAGCCCCCATTGACCAGCTCAAACCCACTAACCATAGAGCTTTCATCGCCCTCTACGTGCAAGTCACACTCATCCTACATTAGAAGCACCGCCCACCATGTACAACTCAGAGCGTCGTAACAAAAAAGGATAACCCTCATGCAGATGCTTTCTTTTCCTACGCTAATCGATCAAATCCTAGATGGGCACCAAACCAGAAACAATGAGACCCACTCTATAACAACCTGTCCCTTATGCCCTCCCATTGCATTCCCCCACGATTACAATTCATTAACGTCCGACGTCTTGCCAGCAGCTGGTAACGAGTCAGCAGCATGCATAACTCCCGTGCTTCTACCACATAGCAGCACGGACTTTATTTCCCCGGTCCACTTACAACTTCATATACTAACTGATTGCTGATCGAGACGTTCCACACTTAATACCAACAAACAGACTTCACATGCGTTTAACCTCACCCCATAAGTATACTTCTAGTACAACGGGCGTCCGCATACTTGCATTCAATTTTTTATCTAGTCAACCATCCTCTATGACTCCCGAATCCAAACTCAAGCACACCCAAAACAGTAAAAAGATTCCGAACACAACTTAATGACCTAACCATAAGTAGGACGCCTCCGTATAAAAACCCTGAGTGAAATCACTGCTATATAACCCTCTCACCTCGTCTTATTACATAAACTAAAATTTCCTACGTTAGTAACGATATCCTGTGCTAGAATACTTGATGAACCCCACTTACCACGCAGCATGTAATCGTTGAATCAAAAGGGGCATTCCTAAAGAACGTTGCTGACAACTGCGATGACACTAATTTGGAAGTCACTCTCACAAAGATCCTCCCCCTGCGCCAGTCCCCCTGTACACATCGCACCCAACCATCCAGAGCTTATAAACAATAAAGTCCTAGTGGAGCCATTCGGCGTGTTCCCACGGTACTATATGGACAGGTCACCACTAATACCTACATGGGTAACACATGCAATCCGCGAGCGTGACACATTTCATCCACCGGTCTCGTTATCAACTTGCCGAGCGTCTTAAGACCAACTAGAAAATCACTATCAACCCTCAGTAGAACGTAAATCCAACTCTTATCAAAACCGTTTGGGGATAACAAACCAAACTACCCAAACGGCGGCAAGCTTTATCTTCAAGAAATGGTATAAGGCTCCTACCATTTCTGCCAATAGACATAGAACTTCTATACATTCCCTAGTCTATTATAATCTAAACTCAAAAGATGGTATCCATAACTTAGAGAGCGAATTCCCACTCATGAAAAGCCTTTCTGCCATTAATCTCTTCAGATAAACCTAACTACGATTATCTCTATTGAAAATGTTCTTACAATCTATTAACTTTCACACAACATCCACTCCTCACCACATTACATCACACCAAACCGAAATCATCTCTATCAAGGTAGCCAGATCTGACTCATCGAAATTACCACACCACTACGTAACATACAAAGACATTCTATATTAACACACTAGAAACAAGCGACAACCACCAACGCAACTGTACCACGGCCAGCCAGCTATTATCCGTAAGAAAAATCATAGTGCGCCGAGCGTTTAACTATAATGCCCTCGCACCATTTATGCGTACACTCTCCGTTTATGTCTCACCGCTGACTCAGTATTATACCATCCAAAGCTGTCTCATTCCCCCAACCCGATTAACTACCGCGCTTACCTATGCAACTACAAATCCAATTCATATAAAACTGGTTACTTCCATGACACGAGATCCACGTCCTACTACTCACCTTTAGTCTTACCGCCTTTACTCACCCCAGTACATATCATGTCCGTACGAGCATCAAACGCCCCTAGAGATAAATACCAAACTAGTATTAGGTAAAATAACCGCACATCTCTAGAATTTCTTAAATTAACATTCCTAAACTCAATAAACGTCGGAAGTAAAACTAGGGTTTGGGTCATATTGTTAGATCAATAATGAGATCTGATACCAAAACATCTAATAACAAAGCCCGTTTAATTCTCAAATACTATCCCCCCTAATTCTCCCCCCATTCCTTCTCTCGTTCTCCTACCAACGAATGCTACGATACCGAACAACGATCTCCCCCAAGATTTCATTAGGCCTCCTCAAATTGTTAAGATTAGCATCAAATACTCGACATCAACCACTCTTGTGAATACCATAACCACTATGCATCGATCCAATAGATCTCAGACTCATTACATAAACTTTCATAATACCCGACTCAACTCTGGATTACCCGACAATCCCACACAAAAGTACCTTTTCTTCCATTCCAACTCTCCCCCCTCTAAACAAATCTCCGAGACAAACACTAGTAGGTACTTCACTCTTCGTATAATATCACAGGTGAGAATCTCACCGCGAACCCTTCTGTCTAATACAGACCCTACAGCACCTTTTCAAATCGTCACACAACCCCGGCCTTTCCCACAAGCCGGATTACCCCAATTAGGATCTATCACGATTATCCACTAGATACCAGCGTACACCTTCGAAAGCTCAGGTAACGAACAGACTTGAATTTAACACTGATAAATAAAACCCCTGCGCCTACCAGTTAGGCTTACCATATGGTCGTAAAGTCTAGTCCCAAACTAACCCTTTCCTAGATCTCCTTCCACACCCCCGCCCCACCGCCCTTCAGACGCGGAATGATACTACAAGTTAATCCTAGTATTCGTGCCTCCCGGTTCACGTACATCTCCCTCCTTTGTATAATTCAGTATTTGGAACCTCACGAAGAGATAGTACTCACGACTCCTAAATACGAGCTTAATAGTACACGATCCTCATAGACTCGAACTGTAAAAACTTTAATGATCACTAGCCCATCATAAACAAAAGGCGGCAATTGGCTCCCCTCGCATACATGACTACCTACTCGTATGCCCACAAACTTTAGCCCACGTACCCTATCCCTCAGCGTATAACTGTCTCATCAAACCAACTCAGAAAAATCTACATTTTCACTCGCTTTAGTAGCCTAAGCTCGAAGCTACTAAGAATTCTCACCTTTTACACATGAAAAGCTAACGGAATTCTTAAAAGACCCGAAGAAAAAGCCAATAAACCAAAAAAAGAGGATCGACGCGTAATCTATCAGCATCTCTTCAACAACCACTTACCATTTCTTTTGCGATAGGAACCCCGCGCACCCCCTAGGAATAACAATCACGTTTAACGCGGATACTACATTTAAACGATCCAACCAGCCATATTTCACACATCGACAATAACCCAACTTAGGCCTACCCACGGAAAACACACATGACTCGATAACATCTCGAACCCATCGCAGCCGGTGTATATGCCCTGTCTCGGCAGTCTGCCCTCATCTACAACTGTACAACATGCTCACAGAGTTTCATCAAATGAAATCGCCAGGATGAAAGGTTCTTCAACCTATGCTTACAGATCACGTACTGCATTAAGCTCGGAAATTAAGGAGCTCATAAGACTACTTCTACGCCACCCACCAATTTATAATATAACAAGAGAACAGACACGCACCCTCATTACACCTCGTTCCCAATCGAGAGTTTCCTACTGTTACCACAGAGTTATTGCATTTCTCACTCTTAACACTATCAATCCATACCCACCAGATACATACGATACCACCTCCTGCCTTTCAGGGCCTACAGATCCCACACCTTCACTTAACATGACAATCACTTAAACAGAAGTCTTACTGGTTGGCGAGATCGTCGACTTCATTGCCACTAAGAATCGTAATTTGAACGATCAATAGGAACGAAAAAGCCTCCCATATTCCGTGTTAATCAACAACAAATCTAACTCGATATCTTAGCAAGCAAAGGGACACCCAGAGCCTTTCCTCCTGGACCACCCCAACCCATGCAATTTTGAAAAAAACTTCCGCGCCGAATTACCCGTTTGCCATAATAATTTAACCAGGCTCTTAGAGAATAACGATTACCAGTCTATCCTCACGTAACGCAAGTCGTGGTACTATTAAACATCCCAAAGCGGCTCCCTGAAGTTCTTTCGCGGCTTTTGCGTATCCACCGCTACGTCTCCCACTAGTCGCAGTACCTTTTATGGCGACAAAAGGATCACACTTTATAGCTGGAGATCGCCTTACATATAACGAAATTTCAAAGGATAAAGAAGATCCCTCACTCGACTCAAAACCAACTTCTTGATAATAGGAAACCGTTTCACCTTACAAGGTATCTCCAAACATGCGACGTCCTTCTCAAGCCGTAATCCCTAGTACCCAGACCAATACCTAATACTACTGTTATAAGATACTTAACTGTCTCTTACAAACTTAAACATTACACACCCCACCCGAGTTATACCGAGTAGCAAAACTGCTGTCAACCACTCAGTATGCACTCCAACTAGCTTAACATAAGTCCCTGAAGCGACAGAAGATCCTACACCAAGTGTCACCTTATCGACGAGGAAAAGCCCGAATTCATGCTGCCAAACATCATAACATAATAACACGATCAACAAAGCCTGATCGCATCTCAGTTAATAATACAACATACAATTACAATGAAAGCCAGAGGCAAAGTCTAGGTATGGCTGAATCCCGCTGAGATAGTCCATCCAAATGCCACCGAGAGATACCACTCAAAATCCCCACCGTCCACCTTCACTCGACACGTTTCATCGCAGGATTAACTTTGTTTTCTACCTTCCATTTTCCATCCCCGTCACGACCCTTCTAGTTCCAACCACACCACTTGTTAGAACATCCTTCACTCTTCTAATAAAAAGCCCTTAAGTATTCTCCGCAGATATGAAACCTCAAGGTATCTATACAGCTAAAATACAGCTGTAAATACAAGAACCAGGTTAACGCAGCCATCACTTCAGCTGACAGGAAGTAACCTTTTCAAAACACGAATCGCCCCACTAGCATGAACAAGTCTAGTTGAGCAGGTTTACAAAGCGTCTCATCATTAAACCAAGCTACTCCCGACCTAAGTCCAGGCACACTTCTAAATCAACAAAGTAAAAAACCAGTGTCTAACGGAAGTCCCAGACCAATCCACACATCACTCAACGCATGCCTCAGCTAAAAGGGACAACCAAGCCGGGCCCACTACGCTCAAATGTAAAAACCACAATGAGTGGCCAAGGGGTATCACACCCAATCCTAATACATACTATGTAGCTTTCCTCCATCGCCAGGTGAATGCTCTGCTTTGATACCTGAAAAATAGCCTTCCGCACGACTATAATTAATAACCTCCATGGACATTCCAGTTAAGAAATTAAACGAATTGAACATGTGAAAAATCAGAGCGATGTCCACCCAAGATTCCTCCGGTCCCTGTTCATTTAATTATTACTAGTTGCAGCTGGGATGTAAAATGCCTCTATCTTAGATCTGACCGAACACGAGCTCGTTTAATCAAAAGTCACTAAAGCGCCCAAACTACCGAATGGCTTTTCCCTGCCTCCACAGAGAAATAAACGTCCCTGTTCATCCTTTCAACCCGCCACCTCAACAATTTAGGATACTGACCTCCGACCCAGTGCGCTTTTAACGAGACTATCATTCTCATGTGTAATCTTATAAAATGCGACCAACTCATAACCATTCCTTCCTCAACTCGTGCTTCCTATACTGGACTACTCCTAGCTTTGGACCTAACAAGAGGCGCTCCGGGCCATTACCACAAGACAAACATTTCGTAAATGAAATAATCTAAAGTGAACTGTTAATGAAGACATTGCTAGCCCATGGTTCCGAGATGGTATGAATTGACTCCAACAACATCATTCCGCACAAACATTTCTAGATCGATGCTCTCATAGGTACCAGGGACCCTCCGAAACAAGTGGGCTATCTTGCGGCGCATGCATTTCCAACGGAACCAACGACCTAGCATATGCACCTCTTCGCAACATTTGCACGCGAAGTCTGACCACCGTCACATCTGAACAAAACCAATTTAAACTACTTAGAAGAAACCCCGGCAAAACCACAATCCCGACGAATAAAATTCGGATAACAACCGAAAACGGCGTACTTAATTCCCGTACAAACTATACAAGACCAACAGATGCAAATATATTGTCACTTTGACTCAAAACCACCGAGTATTCTCCCGACACTGCTCACACTCCGGGATTCTATCATAGAACAGGACCCCTACATAGAAATAAGCTCTCCTAATTAAGCATTAAATAATGGCATCCAAGTACTCTCGCAAATTTAAAGGTGCTATTACATATAGGGACGGCAGACAAAATAAATTTGGTCAAAGTAT

orfs_random_valid = seqshoworfs(hippo_random, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'MinimumLength', 226);
hippo_random = hippo(randperm(length(hippo)))

% Q9 Translate that ORF into an amino acid sequence, using Matlab (Paste the command you used
% into the report, and the first 50 letters of the resulting AA sequence ? not all of it). Paste all of
% the AA sequence into BLAST and report the results including the 5 closest organisms (common
% and Taxonomic name), the name of the protein and its identity score. Create and show a
% multiple alignment of all 5 proteins retrieved via BLAST plus your original animal?s and Briefly
% explain/discuss the findings. 

hippoPORF = nt2aa(hippo(orfs_valid(1).Start(2) : orfs_valid(1).Stop(2)));
hippoPORF

hippoPORF =

MNENLFASFITPTILGLPLVTLIIMFPSILFPAPTRLITNRLVSIQQ*LIQLVSKQIMNIHNHKGQT*TLILISLILFIGSTNLLGLLPHSFTPTTQLSINLGIAIPL*AGTVIIGFRNKTKISLAHFLPQGTPTPLIPMLVIIETISLFIQPIALAVRLTANITAGHLLMHLIGGATLALINISITTALITFIILVLLTALEFAVAIIQAYVFTLLVSLYLHDNT

length(hippoPORF)

ans =

   226

hippoPORF(1:50)

ans =

MNENLFASFITPTILGLPLVTLIIMFPSILFPAPTRLITNRLVSIQQ*LI


% Q10 Retrieve from Genbank a complete sequence of human mtDNA and report its accession number
% and full name. Check if the above protein is found also in the human mtDNA. Briefly motivate
% your conclusions. 

human = fastaread('~/Downloads/homoSapiens.fasta')

human = 

      Header: 'KY408152.1 Homo sapiens isolate 684 mitochondrion, complete genome'
    Sequence: 'GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGAGCCG?'

human = human.Sequence

human =

GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCGCTTTCCACACAGACATCATAACAAAAAATTTCCACCANNNNNNNNNNNNNNNNNNTCTGGCCACAGCACTTAAACACATCTCTGCCAAACCCCAAAAACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTTTATCTTTTGGCGGTATGCACTTTTAACAGTCACCCCCCAACTAACACATTATTTTCCCCTCCCACTCCCATACTACTAATCTCATCAATACAACCCTCGCCCATCCTACCCAGCACACACACACACACCGCTGCTAACCCCATACCCCGAACCAACCAAACCCCAAAGACACCCCCCACAGTTTATGTAGCTTACCTCCTCAAAGCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCCATAAACAAATAGGTTTGGTCCTAGCCTTTCTATTAGCTCTTAGTAAGATTACACATGCAAGCATCCCCGTTCCAGTGAGTTTACCCTCTAAATCACCACGATCAAAAGGGACAAGCATCAAGCACGCAGCAATGCAGCTCAAAACGCTTAGCCTAGCCACACCCCCACGGGAAACAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGGGATTAGATACCCCACTATGCTTAGCCCTAAACCTCAACAGTTAAATCAACAAAACTGCTCGCCAGAACACTACGAGCCACAGCTTAAAACTCAAAGGACCTGGCGGTGCTTCATACCCCTCTAGAGGAGCCTGTTCTGTAATCGATAAACCCCGATCAACCTCACCACCTCTTGCTCAGCCTATATACCGCCATCTTCAGCAAACCCTGATGAAGGCTACAAAGTAAGCGCAAGTACCCACGTAAAGACGTTAGGTCAAGGTGTAGCCCATGAGGTGGCAAGAAATGGGCTACATTTTCTACCCCAGAAAACTACGATAGCCCTTATGAAACTTAAGGGTCGAAGGTGGATTTAGCAGTAAACTGAGAGTAGAGTGCTTAGTTGAACAGGGCCCTGAAGCGCGTACACACCGCCCGTCACCCTCCTCAAGTATACTTCAAAGGACATTTAACTAAAACCCCTACGCATTTATATAGAGGAGACAAGTCGTAACATGGTAAGTGTACTGGAAAGTGCACTTGGACGAACCAGAGTGTAGCTTAACACAAAGCACCCAACTTACACTTAGGAGATTTCAACTTAACTTGACCGCTCTGAGCTAAACCTAGCCCCAAACCCACTCCACCTTACTACCAGACAACCTTAGCCAAACCATTTACCCAAATAAAGTATAGGCGATAGAAATTGAAACCTGGCGCAATAGATATAGTACCGCAAGGGAAAGATGAAAAATTATAGCCAAGCATAATATAGCAAGGACTAACCCCTATACCTTCTGCATAATGAATTAACTAGAAATAACTTTGCAAGGAGAGCCAAAGCTAAGACCCCCGAAACCAGACGAGCTACCTAAGAACAGCTAAAAGAGCACACCCGTCTATGTAGCAAAATAGTGGGAAGATTTATAGGTAGAGGCGACAAACCTACCGAGCCTGGTGATAGCTGGTTGTCCAAGATAGAATCTTAGTTCAACTTTAAATTTGCCCACAGAACCCTCTAAATCCCCTTGTAAATTTAACTGTTAGTCCAAAGAGGAACAGCTCTTTGGACACTAGGAAAAAACCTTGTAGAGAGAGTAAAAAATTTAACACCCATAGTAGGCCTAAAAGCAGCCACCAATTAAGAAAGCGTTCAAGCTCAACACCCACTACCTAAAAAATCCCAAACATATAACTGAACTCCTCACACCCAATTGGACCAATCTATCACCCTATAGAAGAACTAATGTTAGTATAAGTAACATGAAAACATTCTCCTCCGCATAAGCCTGCGTCAGATTAAAACACTGAACTGACAATTAACAGCCCAATATCTACAATCAACCAACAAGTCATTATTACCCTCACTGTCAACCCAACACAGGCATGCTCATAAGGAAAGGTTAAAAAAAGTAAAAGGAACTCGGCAAATCTTACCCCGCCTGTTTACCAAAAACATCACCTCTAGCATCACCAGTATTAGAGGCACCGCCTGCCCAGTGACACATGTTTAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCATAATCACTTGTTCCTTAAATAGGGACCTGTATGAATGGCTCCACGAGGGTTCAGCTGTCTCTTACTTTTAACCAGTGAAATTGACCTGCCCGTGAAGAGGCGGGCATGACACAGCAAGACGAGAAGACCCTATGGAGCTTTAATTTATTAATGCAAACAGTACCTAACAAACCCACAGGTCCTAAACTACCAAACCTGCATTAAAAATTTCGGTTGGGGCGACCTCGGAGCAGAACCCAACCTCCGAGCAGTACATGCTAAGACTTCACCAGTCAAAGCGAACTACTATACTCAATTGATCCAATAACTTGACCAACGGAACAAGTTACCCTAGGGATAACAGCGCAATCCTATTCTAGAGTCCATATCAACAATAGGGTTTACGACCTCGATGTTGGATCAGGACATCCCGATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAAAGTCCTACGTGATCTGAGTTCAGACCGGAGTAATCCAGGTCGGTTTCTATCTACTTCAAATTCCTCCCTGTACGAAAGGACAAGAGAAATAAGGCCTACTTCACAAAGCGCCTTCCCCCGTAAATGATATCATCTCAACTTAGTATTATACCCACACCCACCCAAGAACAGGGTTTGTTAAGATGGCAGAGCCCGGTAATCGCATAAAACTTAAAACTTTACAGTCAGAGGTTCAATTCCTCTTCTTAACAACATACCCATGGCCAACCTCCTACTCCTCATTGTACCCATTCTAATCGCAATGGCATTCCTAATGCTTACCGAACGAAAAATTCTAGGCTATATACAACTACGCAAAGGCCCCAACGTTGTAGGCCCCTACGGGCTACTACAACCCTTCGCTGACGCCATAAAACTCTTCACCAAGGAGCCCCTAAAACCCGCCACATCTACCATCACCCTCTACATCACCGCCCCGACCTTAGCTCTCACCATCGCTCTTCTACTATGAACCCCCCTCCCCATACCCAACCCCCTGGTCAACCTCAACCTAGGCCTCCTATTTATTCTAGCCACCTCTAGCCTAGCCGTTTACTCAATCCTCTGATCAGGGTGAGCATCAAACTCAAACTACGCCCTGATCGGCGCACTGCGAGCAGTAGCCCAAACAATCTCATATGAAGTCACCCTAGCCATCATTCTACTATCAACATTACTAATAAGTGGCTCCTTTAACCTCTCCACCCTTATCACAACACAAGAACACCTCTGATTACTCCTGCCATCATGACCCTTGGCCATAATATGATTTATCTCCACACTAGCAGAGACCAACCGAACCCCCTTCGACCTTGCCGAAGGGGAGTCCGAACTAGTCTCAGGCTTCAACATCGAATACGCCGCAGGCCCCTTCGCCCTATTCTTCATAGCCGAATACACAAACATTATTATAATAAACACCCTCACCACTACAATCTTCCTAGGAACAACATATGACGCACTCTCCCCTGAACTCTACACAACATATTTTGTCACCAAGACCCTACTTCTAACCTCCCTGTTCTTATGAATTCGAACAGCATACCCCCGATTCCGCTACGACCAACTCATACACCTCCTATGAAAAAACTTCCTACCACTCACCCTAGCATTACTTATATGATATGTCTCCATACCCATTACAATCTCCAGCATTCCCCCTCAAACCTAAGAAATATGTCTGATAAAAGAGTTACTTTGATGGAGTAAATAATAGGAGCTTAAACCCCCTTATTTCTAGGACTATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAATTAATCCCCTGGCCCAACCCGTCATCTACTCTACCATCTTTGCAGGCACACTCATCACAGCGCTAAGCTCGCACTGATTTTTTACCTGAGTAGGCCTAGAAATAAACATGCTAGCTTTTATTCCAGTTCTAACCAAAAAAATAAACCCTCGTTCCACAGAAGCTGCCATCAAGTATTTCCTCACGCAAGCAACCGCATCCATAATCCTTCTAATAGCTATCCTCTTCAACAATATACTCTCCGGACAATGAACCATAACCAATACTACCAATCAATACTCATCATTAATAATCATAATGGCTATAGCAATAAAACTAGGAATAGCCCCCTTTCACTTCTGAGTCCCAGAGGTTACCCAAGGCACCCCTCTGACATCCGGCCTGCTTCTTCTCACATGACAAAAACTAGCCCCCATCTCAATCATATACCAAATCTCTCCCTCACTAAACGTAAGCCTTCTCCTCACTCTCTCAATCTTATCCATCATAGCAGGCAGTTGAGGTGGATTAAACCAAACCCAGCTACGCAAAATCTTAGCATACTCCTCAATTACCCACATAGGATGAATAATAGCAGTTCTACCGTACAACCCTAACATAACCATTCTTAATTTAACTATTTATATTATCCTAACTACTACCGCATTCCTACTACTCAACTTAAACTCCAGCACCACGACCCTACTACTATCTCGCACCTGAAACAAGCTAACATGACTAACACCCTTAATTCCATCCACCCTCCTCTCCCTAGGAGGCCTGCCCCCGCTAACCGGCTTTTTGCCCAAATGGGCCATTATCGAAGAATTCACAAAAAACAATAGCCTCATCATCCCCACCATCATAGCCACCATCACCCTCCTTAACCTCTACTTCTACCTACGCCTAATCTACTCCACCTCAATCACACTACTCCCCATATCTAACAACGTAAAAATAAAATGACAGTTTGAACATACAAAACCCACCCCATTCCTCCCCACACTCATCGCCCTTACCACGCTACTCCTACCTATCTCCCCTTTTATACTAATAATCTTATAGAAATTTAGGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTGCAATACTTAATTTCTGTAACAGCTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTACTAGACCAATGGGACTTAAACCCACAAACACTTAGTTAACAGCTAAGCACCCTAATCAACTGGCTTCAATCTACTTCTCCCGCCGCCGGGAAAAAAGGCGGGAGAAGCCCCGGCAGGTTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAATCACCTCGGAGCTGGTAAAAAGAGGCCTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCACTCAGCCATTTTACCTCACCCCCACTGATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGTTTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTATAGTGGAAGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATTTCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACCTATTCTGATTTTTCGGTCACCCTGAAGTTTATATTCTTATCCTACCAGGCTTCGGAATAATCTCCCATATTGTAACTTACTACTCCGGAAAAAAAGAACCATTTGGATACATAGGTATGGTCTGAGCTATGATATCAATTGGCTTCCTAGGGTTTATCGTGTGAGCACACCATATATTTACAGTAGGAATAGACGTAGACACACGAGCATATTTCACCTCCGCTACCATAATCATCGCTATCCCCACCGGCGTCAAAGTATTTAGCTGACTCGCCACACTCCACGGAAGCAATATGAAATGATCTGCTGCAGTGCTCTGAGCCCTAGGATTCATCTTTCTTTTCACCGTAGGTGGCCTGACTGGCATTGTATTAGCAAACTCATCACTAGACATCGTACTACACGACACGTACTACGTTGTAGCTCACTTCCACTATGTCCTATCAATAGGAGCTGTATTTGCCATCATAGGAGGCTTCATTCACTGATTTCCCCTATTCTCAGGCTACACCCTAGACCAAACCTACGCCAAAATCCATTTCACTATCATATTCATCGGCGTAAATCTAACTTTCTTCCCACAACACTTTCTCGGCCTATCCGGAATGCCCCGACGTTACTCGGACTACCCCGATGCATACACCACATGAAACATCCTATCATCTGTAGGCTCATTCATTTCTCTAACAGCAGTAATATTAATAATTTTCATGATTTGAGAAGCCTTCGCTTCGAAGCGAAAAGTCCTAATAGTAGAAGAACCCTCCATAAACCTGGAGTGACTATATGGATGCCCCCCACCCTACCACACATTCGAAGAACCCGTATACATAAAATCTAGACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGGCCTCCATGACTTTTTCAAAAAGGTATTAGAAAAACCATTTCATAACTTTGTCAAAGTTAAATTATAGGCTAAATCCTATATATCTTAATGGCACATGCAGCGCAAGTAGGTCTACAAGACGCTACTTCCCCTATCATAGAAGAGCTTATCACCTTTCATGATCACGCCCTCATAATCATTTTCCTTATCTGCTTCCTAGTCCTGTATGCCCTTTTCCTAACACTCACAACAAAACTAACTAATACTAACATCTCAGACGCTCAGGAAATAGAAACCGTCTGAACTATCCTGCCCGCCATCATCCTAGTCCTCATCGCCCTCCCATCCCTACGCATCCTTTACATAACAGACGAGGTCAACGATCCCTCCCTTACCATCAAATCAATTGGCCACCAATGGTACTGAACCTACGAGTACACCGACTACGGCGGACTAATCTTCAACTCCTACATACTTCCCCCATTATTCCTAGAACCAGGCGACCTGCGACTCCTTGACGTTGACAATCGAGTAGTACTCCCGATTGAAGCCCCCATTCGTATAATAATTACATCACAAGACGTCTTGCACTCATGAGCTGTCCCCACATTAGGCTTAAAAACAGATGCAATTCCCGGACGTCTAAACCAAACCACTTTCACCGCTACACGACCGGGGGTATACTACGGTCAATGCTCTGAAATCTGTGGAGCAAACCACAGTTTCATGCCCATCGTCCTAGAATTAATTCCCCTAAAAATCTTTGAAATAGGGCCCGTATTTACCCTATAGCACCCCCTCTACCCCCTCTAGAGCCCACTGTAAAGCTAACTTAGCATTAACCTTTTAAGTTAAAGATTAAGAGAACCAACACCTCTTTACAGTGAAATGCCCCAACTAAATACTACCGTATGGCCCACCATAATTACCCCCATACTCCTTACACTATTCCTCATCACCCAACTAAAAATATTAAACACAAACTACCACCTACCTCCCTCACCAAAGCCCATAAAAATAAAAAATTATAACAAACCCTGAGAACCAAAATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCCTAGGCCTACCCGCCGCAGTACTGATCATTCTATTTCCCCCTCTATTGATCCCCACCTCCAAATATCTCATCAACAACCGACTAATCACCACCCAACAATGACTAATCAAACTAACCTCAAAACAAATGATAACCATACACAACACTAAAGGACGAACCTGATCTCTTATACTAGTATCCTTAATCATTTTTATTGCCACAACTAACCTCCTCGGACTCCTGCCTCACTCATTTACACCAACCACCCAACTATCTATAAACCTAGCCATGGCCATCCCCTTATGAGCGGGCGCAGTGATTATAGGCTTTCGCTCTAAGATTAAAAATGCCCTAGCCCACTTCTTACCACAAGGCACACCTACACCCCTTATCCCCATACTAGTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTGGCCGTACGCCTAACCGCTAACATTACTGCAGGCCACCTACTCATGCACCTAATTGGAAGCACCACCCTAGCAATATCAACCATTAACCTTCCCTCTACACTTATCATCTTCACAATTCTAATTCTACTGACTATCCTAGAAATCGCTGTCGCCTTAATCCAAGCCTACGTTTTCACACTTCTAGTAAGCCTCTACCTGCACGACAACACATAATGACCCACCAATCACATGCCTATCATATAGTAAAACCCAGCCCATGACCCCTAACAGGGGCCCTCTCAGCCCTCCTAATGACCTCCGGCCTAGCCATGTGATTTCACTTCCACTCCATAACGCTCCTCATACTAGGCCTACTAACCAACACACTAACCATATACCAATGATGGCGCGATGTAACACGAGAAAGCACATACCAAGGCCACCACACACCACCTGTCCAAAAAGGCCTTCGATACGGGATAATCCTATTTATTACCTCAGAAGTTTTTTTCTTCGCAGGATTTTTCTGAGCCTTTTACCACTCCAGCCTAGCCCCTACCCCCCAATTAGGAGGGCACTGGCCCCCAACAGGCATCACCCCGCTAAATCCCCTAGAAGTCCCACTCCTAAACACATCCGTATTACTCGCATCAGGAGTATCAATCACCTGAGCTCACCATAGTCTAATAGAAAACAACCGAAACCAAATAATTCAAGCACTGCTCATTACAATTTTACTGGGTCTCTATTTTACCCTCCTACAAGCCTCAGAGTACTTCGAGTCTCCCTTCACCATTTCCGACGGCATCTACGGCTCAACATTTTTTGTAGCCACAGGCTTCCACGGACTTCACGTCATTATTGGCTCAACTTTCCTCACTATCTGCTTCATCCGCCAACTAATATTTCACTTTACATCCAAACATCACTTTGGCTTCGAAGCCGCCGCCTGATACTGGCATTTTGTAGATGTGGTTTGACTATTTCTGTATGTCTCCATCTATTGATGAGGGTCTTACTCTTTTAGTATAAATAGTACCGTTAACTTCCAATTAACTAGTTTTGACAACATTCAAAAAAGAGTAATAAACTTCGCCTTAATTTTAATAATCAACACCCTCCTAGCCTTACTACTAATAATTATTACATTTTGACTACCACAACTCAACGGCTACATAGAAAAATCCACCCCTTACGAGTGCGGCTTCGACCCTATATCCCCCGCCCGCGTCCCTTTCTCCATAAAATTCTTCTTAGTAGCTATTACCTTCTTATTATTTGATCTAGAAATTGCCCTCCTTTTACCCCTACCATGAGCCCTACAAACAACTAACCTGCCACTAATAGTTATGTCATCCCTCTTATTAATCATCATCCTAGCCCTAAGTCTGGCCTATGAGTGACTACAAAAAGGATTAGACTGAGCCGAATTGGTATATAGTTTAAACAAAACGAATGATTTCGACTCATTAAATTATGATAATCATATTTACCAAATGCCCCTCATTTACATAAATATTATACTAGCATTTACCATCTCACTTCTAGGAATACTAGTATATCGCTCACACCTCATGTCCTCCCTACTATGCCTAGAAGGAATAATACTATCGCTGTTCATTATAGCTACTCTCATAACCCTCAACACCCACTCCCTCTTAGCCAATATTGTGCCTATTGCCATACTAGTCTTTGCCGCCTGCGAAGCAGCGGTGGGCCTAGCCCTACTAGTCTCAATCTCCAACACATATGGCCTAGACTACGTACATAACCTAAACCTACTCCAATGCTAAAACTAATCGTCCCAACAATTATATTACTACCACTGACATGACTTTCCAAAAAACACATAATTTGAATCAACACAACCACCCACAGCCTAATTATTAGCATCATCCCTCTACTATTTTTTAACCAAATCAACAACAACCTATTTAGCTGTTCCCCAACCTTTTCCTCCGACCCCCTAACAACCCCCCTCCTAATACTAACTACCTGACTCCTACCCCTCACAATCATGGCAAGCCAACGCCACTTATCCAGTGAACCACTATCACGAAAAAAACTCTACCTCTCTATACTAATCTCCCTACAAATCTCCTTAATTATAACATTCACAGCCACAGAACTAATCATATTTTATATCTTCTTCGAAACCACACTTATCCCCACCTTGGCTATCATCACCCGATGAGGCAACCAGCCAGAACGCCTGAACGCAGGCACATACTTCCTATTCTACACCCTAGTAGGCTCCCTTCCCCTACTCATCGCACTAATTTACACTCACAACACCCTAGGCTCACTAAACATTCTACTACTCACCCTCACTGCCCAAGAACTATCAAACTCCTGAGCCAACAACTTAATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTTATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGTACTCTTGAAACTAGGCGGCTATGGCATAATACGCCTCACACTCATTCTCAACCCCCTGACAAAACACATAGCCTACCCCTTCCTTGTACTATCCCTATGAGGCATAATTATAACAAGCTCCATCTGCCTACGACAAACAGACCTAAAATCGCTCATTGCATACTCTTCAATCAGCCACATAGCCCTCGTAGTAACAGCCATTCTCATCCAAACCCCCTGAAGCTTCACCGGCGCAGTCATTCTCATAATCGCCCACGGACTTACATCCTCATTACTATTCTGCCTAGCAAACTCAAACTACGAACGCACTCACAGTCGCATCATAATCCTCTCTCAAGGACTTCAAACTCTACTCCCACTAATAGCTTTTTGATGACTTTTAGCAAGCCTCGCTAACCTCGCCTTACCCCCCACTATTAACCTACTGGGAGAACTCTCTGTGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCTACTTACAGGACTCAACATACTAGTCACAGCCCTATACTCCCTCTACATATTTACCACAACACAATGGGGCTCACTCACCCACCACATTAACAACATAAAACCCTCATTCACACGAGAAAACACCCTCATGTTCATACACCTATCCCCCATTCTCCTCCTATCCCTCAACCCCGACATCATTACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTTACGACCCCTTATTTACCGAGAAAGCTCACAAGAACTGCTAACTCATGCCCCCATGTCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGCCCCAAGAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGCACACTACTATAACCACCCTAACCCTAACTTCCCTAATTCCCCCCATCCTTACCACCCTCGTTAACCCTAACAAAAAAAACTCATACCCCCATTATGTAAAATCCATTGTCGCATCCACCTTTATTATCAGTCTCTTCCCCACAACAATATTCATGTGCCTAGACCAAGAAGTTATTATCTCGAACTGACACTGAGCCACAACCCAAACAACCCAGCTCTCCCTAAGCTTCAAACTAGACTACTTCTCCATAATATTCATCCCTGTAGCATTGTTCGTTACATGGTCCATCATAGAATTCTCACTGTGATATATAAACTCAGACCCAAACATTAATCAGTTCTTCAAATATCTACTCATCTTCCTAATTACCATACTAATCTTAGTTACCGCTAACAACCTATTCCAACTGTTCATCGGCTGAGAGGGCGTAGGAATTATATCCTTCTTGCTCATCAGTTGATGATACGCCCGAGCAGATGCCAACACAGCAGCCATTCAAGCAATCCTATACAACCGTATCGGCGATATCGGTTTCATCCTCGCCTTAGCATGATTTATCCTACACTCCAACTCATGAGACCCACAACAAATAGCCCTTCTAAACGCTAATCCAAGCCTCACCCCACTACTAGGCCTCCTCCTAGCAGCAGCAGGCAAATCAGCCCAATTAGGTCTCCACCCCTGACTCCCCTCAGCCATAGAAGGCCCCACCCCAGTCTCAGCCCTACTCCACTCAAGCACTATAGTTGTAGCAGGAATCTTCTTACTCATCCGCTTCCACCCCCTAGCAGAAAATAGCCCACTAATCCAAACTCTAACACTATGCTTAGGCGCTATCACCACTCTGTTCGCAGCAGTCTGCGCCCTTACACAAAATGACATCAAAAAAATCGTAGCCTTCTCCACTTCAAGTCAACTAGGACTCATAATAGTTACAATCGGCATCAACCAACCACACCTAGCATTCCTGCACATCTGTACCCACGCCTTCTTCAAAGCCATACTATTTATGTGCTCCGGGTCCATCATCCACAACCTTAACAATGAACAAGATATTCGAAAAATAGGAGGACTACTCAAAACCATACCTCTCACTTCAACCTCCCTCACCATTGGCAGCCTAGCATTAGCAGGAATACCTTTCCTCACAGGTTTCTACTCCAAAGACCACATCATCGAAACCGCAAACATATCATACACAAACGCCTGAGCCCTATCTATTACTCTCATCGCTACCTCCCTGACAAGCGCCTATAGCACTCGAATAATTCTTCTCACCCTAACAGGTCAACCTCGCTTCCCCACCCTTACTAACATTAACGAAAATAACCCCACCCTACTAAACCCCATTAAACGCCTGGCAGCCGGAAGCCTATTCGCAGGATTTCTCATCACTAACAACATTTCCCCCGCATCCCCCTTCCAAACAACAATCCCCCTCTACCTAAAACTCACAGCCCTCGCTGTCACTTTCCTAGGACTTCTAACAGCCCTAGACCTCAACTACCTAACCAACAAACTTAAAATAAAATCCCCACTATGCACATTTTATTTCTCCAACATACTCGGATTCTACCCTAGCATCACACACCGCACAATCCCCTATCTAGGCCTTCTTACGAGCCAAAACCTGCCCCTACTCCTCCTAGACCTAACCTGACTAGAAAAGCTATTACCTAAAACAATTTCACAGCACCAAATCTCCACCTCCATCATCACCTCAACCCAAAAAGGCATAATTAAACTTTACTTCCTCTCTTTCTTCTTCCCACTCATCCTAACCCTACTCCTAATCACATAACCTATTCCCCCGAGCAATTTCAATTACAATATATACACCAACAAACAATGTTCAACCAGTAACTACTACTAATCAACGCCCATAATCATACAAAGCCCCCGCACCAATAGGATCCTCCCGAATCAACCCTGACCCCTCTCCTTCATAAATTATTCAGCTTCCTACACTATTAAAGTTTACCACAACCACCACCCCATCATACTCTTTCACCCACAGCACCAATCCTACCTCCATCGCTAACCCCACTAAAACACTCACCAAGACCTCAACCCCTGACCCCCATGCCTCAGGATACTCCTCAATAGCCATCGCTGTAGTATATCCAAAGACAACCATCATTCCCCCTAAATAAATTAAAAAAACTATTAAACCCATATAACCTCCCCCAAAATTCAGAATAATAACACACCCGACCACACCGCTAACAATCAATACTAAACCCCCATAAATAGGAGAAGGCTTAGAAGAAAACCCCACAAACCCCATTACTAAACCCACACTCAACAGAAACAAAGCATACATCATTATTCTCGCACGGACTACAACCACGACCAATGATATGAAAAACCATCGTTGTATTTCAACTACAAGAACACCAATGACCCCAATACGCAAAATTAACCCCCTAATAAAATTAATTAACCACTCACTCATCGACCTCCCCACCCCATCCAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAATCACCACAGGACTATTCCTAGCCATGCACTACTCACCAGACGCCTCAACCGCCTTTTCATCAATCGCCCACATCACTCGAGACGTAAATTATGGCTGAATCATCCGCTACCTTCACGCCAATGGCGCCTCAATATTCTTTATCTGCCTCTTCCTACACATCGGGCGAGGCCTATATTACGGATCATTTCTCTACTCAGAAACCTGAAACATCGGCATTATCCTCCTGCTTGCAACTATAGCAACAGCCTTCATAGGCTATGTCCTCCCGTGAGGCCAAATATCATTCTGAGGGGCCACAGTAATTACAAACTTACTATCCGCCATCCCATACATTGGGACAGACCTAGTTCAATGAATCTGAGGAGGCTACTCAGTAGACAGTCCCACCCTCACACGATTCTTTACCTTTCACTTCATCTTGCCCTTCATTATTGCAGCCCTAGCAGCACTCCACCTCCTATTCTTGCACGAAACGGGATCAAACAACCCCCTAGGAATCACCTCCCATTCCGATAAAATCACCTTCCACCCTTACTACACAATCAAAGACGCCCTCGGCTTACTTCTCTTCCTTCTCTCCTTAATGACATTAACACTATTCTCACCAGACCTCCTAGGCGACCCAGACAATTATACCCTAGCCAACCCCTTAAACACCCCTCCCCACATCAAGCCCGAATGATATTTCCTATTCGCCTACACAATTCTCCGATCCGTCCCTAACAAACTAGGAGGCGTCCTTGCCCTATTACTATCCATCCTCATCCTAGCAATAATCCCCATCCTCCATATATCCAAACAACAAAGCATAATATTTCGCCCACTAAGCCAATCACTTTATTGACTCCTAGCCGCAGACCTCCTCATTCTAACCTGAATCGGAGGACAACCAGTAAGCTACCCTTTTACCATCATTGGACAAGTAGCATCCGTACTATACTTCACAACAATCCTAATCCTAATACCAACTATCTCCCTAATTGAAAACAAAATACTCAAATGGACCTGTCCTTGTAGTATAAACTAATACACCAGTCTTGTAAACCGGAGATGAAAACCTTTTTCCAAGGACAAATCAGAGAAAAAGTCTTTAACTCCACCATTAGCACCCAAAGCTAAGATTCTAATTTAAACTATTCTCTGTTCTTTCATGGGGAAGCAGATTTGGGTACCACCCAAGTATTGACTCACCCATCAACAACCGCTATGTATTTCGTACATTACTGCCAGCCACCATGAATATTGTACGGTACCATAAATACTTGACCACCTGTAGTACATAAAAACCCAATCCACATCAAAACCCCCTCCCCATGCTTACAAGCAAGTACAGCAATCAACCCCCAACTATCACACATCAACTGTAACTCCAAAGCCACCCCTCACCCACTAGGATACCAACAAACCTACCCACCCTTAACAGTACATAGCACATAAAGCCATTTACCGTACATAGCACATTACAGTCAAATCCCTTCTCGTCCCCATGGATGACCCCCCTCAGATAGGGGTCCCTTGACCACCATCCTCCGTGAAATCAATATCCCGCACAAGAGTGCTACTCTCCTCGCTCCGGGCCCATAACACTTGGGGGTAGCTAAAGTGAACTGTATCCGACATCTGGTTCCTACTTCAGGGCCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATG

humanProtein = nt2aa(human)
Error using nt2aa (line 165)
The sequence contains the character(s) N.
NT2AA only supports A,C,G,T, and U.
Set the 'ACGTOnly' option to false to resolve ambiguous or unknown characters.Gaps are supported only when a complete codon
is made up of gaps.
 
humanProtein = nt2aa(human, 'ACGTOnly', false)

humanProtein =

DHRSITLLTTHGSSPCIWYFRLGGVHAIALRDAGAGAPYVAVSVFDSCLILLFIAPTFNITGEHTY*SVLIN*CL*DIIITIECLHSRFPHRHHNKKFPPXXXXXXSGHST*THLCQTPKTKNPNTSLTRFQILSFGGMHF*QSPPN*HIIFPSHSHTTNLINTTLAHPTQHTHTHRC*PHTPNQPNPKDTPHSLCSLPPQSNTLKMFRRAHITP*TNRFGPSLSISS**DYTCKHPRSSEFTL*ITTIKRDKHQARSNAAQNA*PSHTPTGNSSD*PLAINESLTKLY*PQGWSISCQPPRSHD*PKSIEAGVKSVLDHPLPNKAKTHLSCKKLQLTQNRLRKWL*HI*THNS*DPNWD*IPHYA*P*TSTVKSTKLLARTLRATA*NSKDLAVLHTPLEEPVL*SINPDQPHHLLLSLYTAIFSKP**RLQSKRKYPRKDVRSRCSP*GGKKWATFSTPENYDSPYET*GSKVDLAVN*E*SA*LNRALKRVHTARHPPQVYFKGHLTKTPTHLYRGDKS*HGKCTGKCTWTNQSVA*HKAPNLHLGDFNLT*PL*AKPSPKPTPPYYQTTLAKPFTQIKYRR*KLKPGAIDIVPQGKDEKL*PSII*QGLTPIPSA**IN*K*LCKESQS*DPRNQTSYLRTAKRAHPSM*QNSGKIYR*RRQTYRAW**LVVQDRILVQL*ICPQNPLNPLVNLTVSPKRNSSLDTRKKPCRESKKFNTHSRPKSSHQLRKRSSSTPTT*KIPNI*LNSSHPIGPIYHPIEELMLV*VT*KHSPPHKPASD*NTELTINSPISTINQQVIITLTVNPTQACS*GKVKKSKRNSANLTPPVYQKHHL*HHQY*RHRLPSDTCLTAAVP*PCKGSIITCSLNRDLYEWLHEGSAVSYF*PVKLTCP*RGGHDTARREDPMEL*FINANST*QTHRS*TTKPALKISVGATSEQNPTSEQYMLRLHQSKRTTILN*SNNLTNGTSYPRDNSAILF*SPYQQ*GLRPRCWIRTSRWCSRY*RFVCSTIKVLRDLSSDRSNPGRFLSTSNSSLYERTREIRPTSQSAFPRK*YHLNLVLYPHPPKNRVC*DGRAR*SHKT*NFTVRGSIPLLNNIPMANLLLLIVPILIAMAFLMLTERKILGYIQLRKGPNVVGPYGLLQPFADAIKLFTKEPLKPATSTITLYITAPTLALTIALLL*TPLPIPNPLVNLNLGLLFILATSSLAVYSIL*SG*ASNSNYALIGALRAVAQTISYEVTLAIILLSTLLISGSFNLSTLITTQEHL*LLLPS*PLAII*FISTLAETNRTPFDLAEGESELVSGFNIEYAAGPFALFFIAEYTNIIIINTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFL*IRTAYPRFRYDQLIHLL*KNFLPLTLALLI*YVSIPITISSIPPQT*EICLIKELL*WSK**ELKPPYF*DYENRTHP*ESKILRATYHTPS*SKVS*ISYRAHTPKMLVIPFPY*LIPWPNPSSTLPSLQAHSSQR*ARTDFLPE*A*K*TC*LLFQF*PKK*TLVPQKLPSSISSRKQPHP*SF**LSSSTIYSPDNEP*PILPINTHH**S*WL*Q*N*E*PPFTSESQRLPKAPL*HPACFFSHDKN*PPSQSYTKSLPH*T*AFSSLSQSYPS*QAVEVD*TKPSYAKS*HTPQLPT*DE**QFYRTTLT*PFLI*LFILS*LLPHSYYST*TPAPRPYYYLAPETS*HD*HP*FHPPSSP*EACPR*PAFCPNGPLSKNSQKTIASSSPPS*PPSPSLTSTSTYA*STPPQSHYSPYLTT*K*NDSLNIQNPPHSSPHSSPLPRYSYLSPLLY**SYRNLG*IQTKSLQSPQ*VAILNFCNS*GLQNPTLHQLNANQPL*LS*ALTRPMGLKPTNT*LTAKHPNQLASIYFSRRREKRREKPRQV*SCFFEFAIQYENHLGAGKKRPNPCL*IYSPMLHSAILPHPH*CSPTVDYSLQTTKTLEHYTYYSAHELES*AQL*ASLFEPSWASQATF*VTTTSTTLSSQPMHL**SSS**YPS*SEALATD*FP**SVPPIWRFPA*TT*ASDSYLPLSYSCSHLL*WKPEQEQVEQSTLP*QGTTPTLEPP*T*PSSPYT*QVSPLS*GPSISSQQLSI*NPLP*PNTKRPSSSDPS*SQQSYFSYLSQS*LLASLYY*QTATSTPPSSTPPEEETPFYTNTYSDFSVTLKFIFLSYQASE*SPIL*LTTPEKKNHLDT*VWSEL*YQLAS*GLSCEHTIYLQ*E*T*THEHISPPLP*SSLSPPASKYLADSPHSTEAI*NDLLQCSEP*DSSFFSP*VA*LALY*QTHH*TSYYTTRTTL*LTSTMSYQ*ELYLPS*EASFTDFPYSQATP*TKPTPKSISLSYSSA*I*LSSHNTFSAYPECPDVTRTTPMHTPHETSYHL*AHSFL*QQ*Y**FS*FEKPSLRSEKS***KNPP*TWSDYMDAPHPTTHSKNPYT*NLDKKGRNRTPQSWFQANPMASMTFSKRY*KNHFITLSKLNYRLNPIYLNGTCSASRSTRRYFPYHRRAYHLS*SRPHNHFPYLLPSPVCPFPNTHNKTN*Y*HLRRSGNRNRLNYPARHHPSPHRPPIPTHPLHNRRGQRSLPYHQINWPPMVLNLRVHRLRRTNLQLLHTSPIIPRTRRPATP*R*QSSSTPD*SPHSYNNYITRRLALMSCPHIRLKNRCNSRTSKPNHFHRYTTGGILRSML*NLWSKPQFHAHRPRINSPKNL*NRARIYPIAPPLPPLEPTVKLT*H*PFKLKIKRTNTSLQ*NAPTKYYRMAHHNYPHTPYTIPHHPTKNIKHKLPPTSLTKAHKNKKL*QTLRTKMNENLFASFIAPTILGLPAAVLIILFPPLLIPTSKYLINNRLITTQQ*LIKLTSKQMITIHNTKGRT*SLILVSLIIFIATTNLLGLLPHSFTPTTQLSINLAMAIPL*AGAVIIGFRSKIKNALAHFLPQGTPTPLIPILVIIETISLLIQPIALAVRLTANITAGHLLMHLIGSTTLAISTINLPSTLIIFTILILLTILEIAVALIQAYVFTLLVSLYLHDNT**PTNHMPII**NPAHDP*QGPSQPS**PPA*PCDFTSTP*RSSY*AY*PTH*PYTNDGAM*HEKAHTKATTHHLSKKAFDTG*SYLLPQKFFSSQDFSEPFTTPA*PLPPN*EGTGPQQASPR*IP*KSHS*THPYYSHQEYQSPELTIV**KTTETK*FKHCSLQFYWVSILPSYKPQSTSSLPSPFPTASTAQHFL*PQASTDFTSLLAQLSSLSASSAN*YFTLHPNITLASKPPPDTGIL*MWFDYFCMSPSIDEGLTLLV*IVPLTSN*LVLTTFKKE**TSP*F**STPS*PYY**LLHFDYHNSTAT*KNPPLTSAASTLYPPPASLSP*NSS**LLPSYYLI*KLPSFYPYHEPYKQLTCH**LCHPSY*SSS*P*VWPMSDYKKD*TEPNWYIV*TKRMISTH*IMIIIFTKCPSFT*ILY*HLPSHF*EY*YIAHTSCPPYYA*KE*YYRCSL*LLS*PSTPTPS*PILCLLPY*SLPPAKQRWA*PY*SQSPTHMA*TTYIT*TYSNAKTNRPNNYITTTDMTFQKTHNLNQHNHPQPNY*HHPSTIF*PNQQQPI*LFPNLFLRPPNNPPPNTNYLTPTPHNHGKPTPLIQ*TTITKKTLPLYTNLPTNLLNYNIHSHRTNHILYLLRNHTYPHLGYHHPMRQPARTPERRHILPILHPSRLPSPTHRTNLHSQHPRLTKHSTTHPHCPRTIKLLSQQLNMTSLHNSFYSKDTSLRTPLMTP*SPCRSPHRWVNSTCRSTLETRRLWHNTPHTHSQPPDKTHSLPLPCTIPMRHNYNKLHLPTTNRPKIAHCILFNQPHSPRSNSHSHPNPLKLHRRSHSHNRPRTYILITILPSKLKLRTHSQSHHNPLSRTSNSTPTNSFLMTFSKPR*PRLTPHY*PTGRTLCASNHVLLIKYHSPTYRTQHTSHSPILPLHIYHNTMGLTHPPH*QHKTLIHTRKHPHVHTPIPHSPPIPQPRHHYRVFLL*I*FNQNIRL*I*QQRLTTPYLPRKLTRTANSCPHV*QHGFLNF*RITAIHWS*APRILVQLQIKVITMHTTITTLTLTSLIPPILTTLVNPNKKNSYPHYVKSIVASTFIISLFPTTIFMCLDQEVIISN*H*ATTQTTQLSLSFKLDYFSIIFIPVALFVTWSIIEFSL*YINSDPNINQFFKYLLIFLITILILVTANNLFQLFIG*EGVGIISFLLIS**YARADANTAAIQAILYNRIGDIGFILALA*FILHSNS*DPQQIALLNANPSLTPLLGLLLAAAGKSAQLGLHP*LPSAIEGPTPVSALLHSSTIVVAGIFLLIRFHPLAENSPLIQTLTLCLGAITTLFAAVCALTQNDIKKIVAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFMCSGSIIHNLNNEQDIRKIGGLLKTIPLTSTSLTIGSLALAGIPFLTGFYSKDHIIETANISYTNA*ALSITLIATSLTSAYSTRIILLTLTGQPRFPTLTNINENNPTLLNPIKRLAAGSLFAGFLITNNISPASPFQTTIPLYLKLTALAVTFLGLLTALDLNYLTNKLKIKSPLCTFYFSNILGFYPSITHRTIPYLGLLTSQNLPLLLLDLT*LEKLLPKTISQHQISTSIITSTQKGIIKLYFLSFFFPLILTLLLIT*PIPPSNFNYNIYTNKQCSTSNYY*STPIIIQSPRTNRILPNQP*PLSFINYSASYTIKVYHNHHPIILFHPQHQSYLHR*PH*NTHQDLNP*PPCLRILLNSHRCSISKDNHHSP*IN*KNY*THITSPKIQNNNTPDHTANNQY*TPINRRRLRRKPHKPHY*THTQQKQSIHHYSRTDYNHDQ*YEKPSLYFNYKNTNDPNTQN*PPNKIN*PLTHRPPHPIQHLRMMKLRLTPWRLPDPPNHHRTIPSHALLTRRLNRLFINRPHHSRRKLWLNHPLPSRQWRLNILYLPLPTHRARPILRIISLLRNLKHRHYPPACNYSNSLHRLCPPVRPNIILRGHSNYKLTIRHPIHWDRPSSMNLRRLLSRQSHPHTILYLSLHLALHYCSPSSTPPPILARNGIKQPPRNHLPFR*NHLPPLLHNQRRPRLTSLPSLLNDINTILTRPPRRPRQLYPSQPLKHPSPHQARMIFPIRLHNSPIRP*QTRRRPCPITIHPHPSNNPHPPYIQTTKHNISPTKPITLLTPSRRPPHSNLNRRTTSKLPFYHHWTSSIRTILHNNPNPNTNYLPN*KQNTQMDLSL*YKLIHQSCKPEMKTFFQGQIREKVFNSTISTQS*DSNLNYSLFFHGEADLGTTQVLTHPSTTAMYFVHYCQPP*ILYGTINT*PPVVHKNPIHIKTPSPCLQASTAINPQLSHINCNSKATPHPLGYQQTYPPLTVHST*SHLPYIAHYSQIPSRPHG*PPSDRGPLTTILREINIPHKSATLLAPGP*HLGVAKVNCIRHLVPTSGP*SLNSPHVPLK*DITM

orfs_human = seqshoworfs(human, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 226);
humanPORF = nt2aa(human(orfs_human(1).Start(1) : orfs_human(1).Stop(1)));
length(humanPORF)

ans =

   316

humanPORF

humanPORF =

MANLLLLIVPILIAMAFLMLTERKILGYIQLRKGPNVVGPYGLLQPFADAIKLFTKEPLKPATSTITLYITAPTLALTIALLL*TPLPIPNPLVNLNLGLLFILATSSLAVYSIL*SG*ASNSNYALIGALRAVAQTISYEVTLAIILLSTLLISGSFNLSTLITTQEHL*LLLPS*PLAII*FISTLAETNRTPFDLAEGESELVSGFNIEYAAGPFALFFIAEYTNIIIINTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFL*IRTAYPRFRYDQLIHLL*KNFLPLTLALLI*YVSIPITISSIPPQT

[score, globalAlignment] = nwalign(
 [score, globalAlignment] = nwalign(
                                    ?
Error: Expression or statement is incorrect--possibly unbalanced (, {, or [.
 
Did you mean:
hippoProtein = nt2aa(hippo)
hippoProtein =

VNVAQTPKARH*KCLDGLTQPRKQTGLVPAFLLIFNRITHASIHTPVRMPSKSPGSKGAGIKHTIHRSSKRPAQPHPHGKQQ*PKLSHERKFD*AILTRVGKSRASHRGHTIDSN*QKYGVKRVKELKVQIKLNSN*TVKSPN*NKNLLQK*L*QY*PHDS*DPNWD*IPHYA*P*TQIIPKTKLFARVLLATA*NSKDLAVLHTPLEEPVL*SINPDKPHQPLLIQSIYRHLQQTLKRTKSKLNYYT*RR*VKV*PMGWEEMGYIF*NKNTTHPNENSYES*ELKEDLVVNQE*SA*LNKAMKHAHTARHPPQITQTHKHNS*TKPNERRQVVTR*AYWKVCLDKSKHSLDLKHLVYTQKISQ*VNALN*S*PNYPTPL*IKQSI*SPT*SIGDRNPRPFGAIEIVP*GNDERSNKSIKKQRSHLVPFA**FN*QNLSKENLS*ATRNQTSYL*TVYQNKLIYVAK**KDL*VEVKSLTSLVIAGCPEKESKFNFKITKSPNQA*CNFKI*SKKVQLFRQRIQP*LVSKFKQHHSWPKSSHQLR*RSSSTISTT*FPYQANNS*PNYWTILFKHRSNTVNMSNKKQFSQHKPTSATDNILTVNNK*I*PNTKTFTPYTVNPTQACTQGKIKKSKRNSANTNPACLPKTSPLASLVLEALPAQ*QLLNGRGILTVQR*HNHLFLK*GLV*MATRGFYCLLLSISEIDLPVKRRE*YNKTRRPCGASIS*LNKNKINPQGTIKSYMSQ*FWLG*PRRKKNPPSDKNLDSPVKT*QHSLTQNL*STEQVTPGITAQSYSRVHIDNRVYDLDVGSGHPNGAAAIKGSFVQRLKSYVI*VQTGAIQVSFYLLYISPSTKGQEK*GLLPISALKLINDIVLT*LNSINIPALDQGTVAMAEPGNCIKLKPLHQRFKSSSQQNVYHQHPYTCRTHPSSHSIPNTSRTKNPRVYTTPKRPKHCRTIRPTTTLRRCNQALHKRTPTTIHLLCLHIHHRTNPSSNSSPNNMNPPTHTIPPHQHKLGRTIYTSHI*PSCIFHPMIRMSLQLEICINWRPTSSSTNNLI*SNPSNHPPVHSPNKRVLHIINPHHNTGKTVTNLPFMTTGHNMIHLDPSRD*PSPIRPHRRRIRTCIRLQRRVCSRAIRYILHSRIHQHHHNKCLHNSPIPRCIPQPILTRTLHNQLYHQNTTTNNILPVNPSILPTIPIRPTNAPSMKKLPTINTSPVHMTRITPHYNIKHPPSNIRNMSDKRITLIE*IIEVQALLFLEL*ESNLLLRTQNSLCYQSHPVPQ*GQLNKLSGPYPENVGLNPSRTNKSLRLCHYLHNHHPGDHDCNN*LSLTINLNRVRNKHASYYPHHNEITQPTSHRSLR*ILHNPSHRLHATHTSSHYQSTILRTMNSHKNT*SNSIYNYNDSSCHKIRTIPLPLLSP*SNTRHPLNSRLNLTDMTKTRTPIHSLPNLPINQSKPNPNYINAIHPSRRMRRAKPNSITKNHSVLIYCPYGVDSSHPNLQPNHNHPKPNNLPHNNLHNVYNICTQLNHHYPFPITHMKQNPHYYNPYTHYSTINRRTTPTNRLRTKMNNHPRNNKKR*HYPTHTNSHHSTPQPIFLHTTHLLYSTNHISFLKQHKNKMTI*SLKTQNTPANNNHPLHHTPSPHTNTSSTRLEV*VT*TKSLQSSKQVYIT*SLLNRDCKTAPHIN*MQLNCFY*AKSLLDWWGMLPTNF*LTAKYPNQLASIYFSRRGEIKAGEAPAELKLLL*ICNSI*IFHYRTWQKEDSTPVLRFTV*CLLSHFTHVHKPLTILNQPQRHRYTISTIRRLSWHSRHWPEPTNPRRTGSTWHTIRR*PNLQRSCHSPRICNNFLYSYTNYDWRVRKLTCSTNNRSP*YGLSSNK*HKLLTTPSLLPTTISILHGRSRGRNRLDRLSPFSRKFSPCWSLCRSYNFLPPPSRSLLYSRCNQLHYHHHQHETTRYISVSNPTVCLISPNHGCATSTLPTCFSSRYYYATHRSKPKYHLL*PCRRRRPCPLSTPILILRTPRSIHPDPPRLRNNLAHCNILLRKKRTFWVHRHSLSYNIHRVPRIYCMSPSHIYSRYRRRHPSILHIRHYNYRHPHRSKSIQLTSNTAWREHQMVPCYDVSPRLYFPIHSGWPNRYCFSQLIPRHRPSRHLLRSSPFPLRALNRRRLRYHRRLRTLIPTIFRIYTQ*HMSKNPLRNHVRGGQSNFLPTAFLRPIRNAPTILRLPRRLYNMKHYLLNRFFYLTNSCSTNSVHHLRGICLQTRSLGCRSNYNQLRVTKRVPSTIPHI*RTRIRELN*PKQERKESNLLLLVSSQHHNHYVSLHKRGISKNYITSSKLSYR*KPCLPPWHIPSN*AFKMQYHPL*KNYCIFTTTR**SYS*SAH*SFTLLH*Y*LPN*PTQTP*MHKR*KLSEQSYQPLSLS*LHCHLCESSI**TKLTTPP*P*KLWATNDTEVTSIQIMKT*TLTPT*SQHQT*SRGTYDS*K*ITESSYP*M*QFEY*SHQKTYYTPEPCHH*V*KQMPFQDD*TKQP*YQHDPDYFTDSAPKSVVPTTALCPLS*N*SHCKLSKNEPHPYYRLIKKLVALTF*VRD*EPSSP**HATTRHINMIYHHPIHISDPIYYLSTENLKTHLPPKP*DYSSHNTKTAYPLRNEMNENLFASFITPTILGLPLVTLIIMFPSILFPAPTRLITNRLVSIQQ*LIQLVSKQIMNIHNHKGQT*TLILISLILFIGSTNLLGLLPHSFTPTTQLSINLGIAIPL*AGTVIIGFRNKTKISLAHFLPQGTPTPLIPMLVIIETISLFIQPIALAVRLTANITAGHLLMHLIGGATLALINISITTALITFIILVLLTALEFAVAIIQAYVFTLLVSLYLHDNT**PTKPTHTI**TQVPDLLQEPSQPY**RRA*PYDSTLTPLSY*RQD*LPIS*QYISDDEM*SEKAPFKATTHQSYKKDFATE*SYLLSPKSYFSQASSEPFTTQASLLLLN*ADVDHPQASTL*TH*RCHF*TPPFY*PLVSPLPEPTTV**KAIENKYSKPSSSQSP*VCTSHYCKLQNTMKPPLQSQMGFTAQLSL*PQAFMDYM*LLAPLS*LYASYAN*NSTSRQITTLASRPPPDTDTS*M*SDYSSTCPFIDEVHSSFSIKTVQLTSNQLASVHSGKEQLT***HY*QTPH*PLYWSSLPSDSHN*TPTQKKQAPMNADLTL*DQPAYLSL*NSS**PSHSFFST*RSPSYFLSHGQPKQQT*KPYSL*PSP*SHS*QSA*PTNELKRD*NEPNMVFSLKQNK*FRLIKL*TNS*IPSVFSIYKYYHGLYNIPRRTVNISIPLNILTPMSRRNDIITIYHSNSHHPKCTLHPSQHNANYSTSFRSM*SSPRTIATSNGIKHIRYRLRTKPKPSPMLKYIIPTIILIPLT*ISKNSII*TNTTAHSLLISFTSLLLLNQFNDNSLNFSPMFFSDPLSTPLLILTIWLLPLILIASQSHLLKEPPTRKKLFITILVTLQTFLIITFSAIELILFYILFEATLIPTLIIITR*GNQTERLNAGLYFLFYTLIGSLPLLVALIYIQNITGSLNFLILQY*TQAVSNS*SNVFL*LACIIAFMVKIPLYGLHL*LPKAHVEAPIAGSIVLAAILLKLGGYGILRITTILNPLTEIMAYPFIILSL*GIIMTSSICLRQTDLKSLIAYSSVSHIALVIVAILIQTP*SYIGATALIIAHGLTSSILFCLANSNYERIHSRTIILARGLQTLLPLIAA**LLASLTNLALPPSINLVGELLVIMSSFS*SNITIILIGTNIIITALYSLYILTTTQRGKYTHHINNITPSFTRENALIALHILPLLLLSLNPKIILGPLYCKHSLRKTLDCESTNRSSNPSYLPRKHARTANSCPHI*QYGFPRLLKDGSYPLVLGTKKLVQLQIKAINLFSSTTLTILFVLTLPIIITNTNIYKSDKYPTYVKNTVSSAFLISLVPIIAFTNTGQEIIISN*H*ITIQTLKLTLSFKADYFSIVFAPVALFVTWSIMEFSIWYIHSDPHINQFFKYLLLFLITIIVLVTANNLFQLFIG*EGVGIMSFLLIG**HGRTDANTAAIQAILYNRIGDVGFIIAIA*FLSNLNT*DIQQIFIINPTHSNLPLIGLILAATGKSAQFGLHP*LPSAMEGPTPVSALLHSSTIVVAGVFLLIRFYPMIENNNLIQTITICLGAITTLFTAICALTQNDIKKIIAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFICSGSIIHNLNNEQDIRKIGGLFKTIPFTTTTLIVGSIALTGVPFLTGFYSKDLIIEAANTSYSNA*ALLITLVATSLTAVYSTRIIFFALLGHPRFPTSTLINENNPLLLNSLKRLIAGSIFAGFILSHNLPPITTPLITIPPYLKMTALAVTILGFTLAFEITLNTQNLKHKHPTNSFKFSTLLGYFPTIIHRLPPHLSLTASQKLASSLLDSA*LENILPKSIAHAQLKLSTLVSNQKGLIKIYFLSFLITIPLSLILFNPHA*LP*SPQHQQTRTNPLQQPIKHHSYIKLQHP*PPH*KPQNHSYHKQPNHPDH*T*TQAPPLLPSRHIKPLKTLSSTLKETPPTQPY*RLKPRDTAQ*PLQSYNQIPPTSHPNKPKTPLNPKKTRQNSTPYHNQPHH*QSNPAHHKSAKALKKTQ*NQSQRQYSK*IQHTLSLFPYGL*P*PTA*KTIVVIQL*EH**QTSENLTP**KLSTMHSLTSQLHQTSHRDETSAPYLASA*SYKF*QAYSWPYTTHQIHSPHSHR*PTSAVM*TTGESSATYTQTAHPSSSSASLLT*DAAYTMAPTHS*KPETSELSYYSQP*LPRL*ATYCHEDKCHSEGQQSLPTYCQLSPILEQT**NESEEAFP*TKPPLHDSLPSTLFFHSLSQH*PSSIYYSSMKQDPTTQQESPQTQTKSHSTPITQSRTS*VSYS**QHYSH*PYLPQTS*GTQTTTPPQTPLAHHHTLNQNDISCSRTRFSDQSPTN*EAS*P*LSQS*SWP*SQYYTHPNNEA*YFDPSANACFEH*SPTY*HSHELEDNPSNTPSSSSDKSPQSYISS*S*Y*CP*QALSKTNS*NEESL*YMTLPRSCKPKKEATHLPETQGRSSSSTISTQS*NSR*TIP*FPCMYYLQDYKVLT*YYNLKCTYIHGYVRRALLLYHIQ*YILCIIVHRTYYV*SCITLSSTCL*ACTSELLAPHGTC*FFIVHGTCTQIISSHQAYPAP*ITSLITRPRETSNPLGRFPSSRSGPIACGGF*E*TLSGIWFLLQGHLI*NRSFFPLK*DISMD**LISPCSRITVISCIWYFYFLGDAWTQL*PPGLNLVN*LVAGL*LNLIYQPGSELQGLFSQWLQDIEKCDVHSFATGHAFIHCSLALARHAYYSVNGYRTYQYISPPCSFKSLIHLLKYASPKQRAIPLDSTNEFAKN*ISKPAQTTQGSPTS*T

[score, globalAlignment] = nwalign(hippoProtein, humanProtein)

score =

  809.3333


globalAlignment =

----VN-V-AQ--TPKA---R----H*KCL-D-GLTQPR---K--QTG---L-V-PAF-L---------LI---F-N-RIT----HASI----HT----P------------VRM---P-SKSPG-S--K----G-AGI---K-----HTI---H-RSS---K----RPAQ---PH---PH-GKQQ*PK-----L-S---H-ER-K-FD*A-IL--T-RVGKS-R-AS-------H-R-GH------TI--D-----SN-*Q----KY---G-VK----RVKE-L-KV---Q---IKLN---S-N*--T----VKS------PN-----*N-KNL-L-Q----K*L*QY*PHDS*DPNWD*IPHYA*P*TQIIPKTKLFARVLLATA*NSKDLAVLHTPLEEPVL*SINPDKPHQPLLIQSIYRHL-QQTLKRTKSKLNYYT*RR*VKV*PMGWEEMGYIF*NKNTTHPNENSYES*ELKEDLVVNQE*SA*LNKAMKHAHTARHPPQITQTHKHNS*TKPNERRQVVTR*AYWKVCLDK-SKHSLDLK-HLVYTQKISQ*VNAL-N-*S*PN-YPTPL*IKQSI*SPT*SIGDRN-P-RPFGAIEIVP*GNDERSNKSIKKQRSHLVPFA**FN*QNLSKENLS*ATRNQTSYL*TVYQNKLIYVAK**KDL*VEVKSLTSLVIAGCPEKESKFNFKI-TKSPNQA*CNFKI*SKKVQLFRQRIQP*LVS-KFKQHHSWPKSSHQLR*RSSSTISTT*FPYQ-ANNS*PNYWTILFK-HRSNTVNMSNKKQFSQHKPTSATDNILTVNNK*I*PNTKTFTPYTVNPTQACTQGKIKKSKRNSANTNPACLPKTSP-LASLVLEALPAQ*QLLNGRGILTVQR*HNHLFLK*GLV*MATRGFYCLLLSISEIDLPVKRRE*YNKTRRPCGASIS*LNKNKI-NPQGTIKSYMSQ*FWLG*PRRKKNPPSDKNLDSPVKT*QHSLTQNL*STEQVTPGITAQSYSRVHIDNRVYDLDVGSGHPNGAAAIKGSFVQRLKSYVI*VQTGAIQVSFYLLYISPSTKGQEK*GLLPISALKLINDIVLT*LNSINIPALDQGTVAMA-EPGNCIKL-KPLHQRFKSSSQQNVYHQH-PYTCR-THPSSHSIPNTSRTKNPRVYTTPKRPKHCRTIRPTTTLRRCNQALHKRTPTTIHLLCLHIHHRTNPSSNSSPNNMNPPTHTIPPH--QHKLGRTIYTSHI*PSCIFHPMIRMSLQLEICINWRPTSSSTNNLI*SN-PSN-HPPVHSPNKRVLHIINPHHNTGKTVTNLPFMTTG-HNMIHLDPSRD*PSPIRPHRRRIRTCIRLQRRVCSRAIRYILHSRIHQHHHNKCLHNSPIPRCI-PQPILTRTLHNQLYHQNTTTNNILPVNPSILPTIPIRPTNAPSMKKLPTINTSPVHMTRITPHYNIKHPPSNI-RNMSDKRITLIE*IIEVQALL-FL--EL*ESNLLLRTQ-NSLCYQSHPVPQ-*GQLNKLSGPYPENVGLNPSRTNKS-LRLC-HYLHNHHPGDH--DCNN*LSLTIN-LNRVRN-K-HASYYPHHNEITQPTSHRSLR*ILHNPSHRLHATHTSSHYQSTILRTMNSHK-NT-*SNSIYNYNDSSCHKIRTIPLPLLSP*SNTRHPLNSR-LNLTDMT-K-TRTPI-HSLPNLPINQS-KPNPNYINAI-HPSRRMRRAKPNSITKNHSVLIY-C-PYGVDSSH-PNLQP-NH-NHPKPNNLPHNNLHNVYNICTQLNHHYPFPITHMKQNPHYY-NPYT--H--YSTINRRTTPTNRL-R-TKM-NNHPRNNKKR*HYPTHTNSHHS-TPQPIFLHTTHLLYSTNHISFLKQHKNKMT-I*SLKTQNTPANNNH--PL-HHTP-SPHTNTSSTRLEV*VT*TKSLQSSKQVYIT*SLLNRDCKTAPHIN*MQLNCFY*AK-SLLDWWGMLPTNF*LTAKYPNQLASIYFSRRGEIKAGEAPAELKLLL*ICNSI*IFHYRTWQKEDSTPVLRFTV*CLLSHFTHVHKPLTILNQPQRHRYTISTIRRLSWHSRHWPEPTNPRRTGSTWH-TIRR*PNLQRSCHSPRICNNFLYSYTNYDWRVRKLTCS-TNNRSP*YGLSSN-K*HKLLT--TPSLLPTTISILHGRSRGRNRLDRLSPFSRKFSPCWSLCRSYNFLPPPS-RSLLYS-RCNQLHYHHHQHETTRYISVSNPTVCLISPNHGCATSTLPTCFSSRYYYATHRSKPKYHL-L*PCRRRRPCPLSTPILILRTP-RS-IHPD-PPRLR-NNLAHCNILLRK-K-RTFWVHRHSLSYNI-HRV-PRIYCMSPSHIYSRYRRRHPSILHIRHYNYRHPHRSKSIQLTSNTAWREHQMVPC-YDVSPRLYFPIHSG-W-PNRYCFSQLIPRHR-PSRHLLRSSPFPLRALNRR-RLRYHRRLRTLIPTIFRIYTQ*-HMS-KN-PL-RNHV-RGGQSNFLPTAFLR-PIRNAPTILRLPRRL-YNM-KHYLLNRFFYLTNS-CSTNSVHHLRGICLQTRSLGCRSNYNQLRVTKRVPSTIPHI*RTRIRELN*PKQERKESN-LL-LLVSSQH-HNHYVSLHKRGISKN--YITSSKLSYR*KPCLPPWHIPSN*AFKMQ-YHPL*KNYCIFTTTR**SYS*SAH*SFTLLH*Y*LPN*PTQTP*MH-KR*KLSEQSYQPLSLS*LHCHLCESSI**TKLTTPP*P*KLWATNDTEVTSIQIMKT*TLTPT*SQHQT*SRGTYDS*K*ITESSYP*M*QFEY*SHQKTY-YTPEPCHH*V*KQMPFQDD*TKQP*YQHDPDYFTDSAPKSVVPTTALCPLS*N*SHCKLSK-NEPHPYYRLIKKLV-ALTF*VR-D*EPSSP**HATTRHINMIYH-HPIHISDPIYYLSTENLKTHLPPKP*DYSSHNTKTAYPLRNEMNENLFASFITPTILGLPLVTLIIMFPSILFPAPTRLITNRLVSIQQ*LIQLVSKQIMNIHNHKGQT*TLILISLILFIGSTNLLGLLPHSFTPTTQLSINLGIAIPL*AGTVIIGFRNKTKISLAHFLPQGTPTPLIPMLVIIETISLFIQPIALAVRLTANITAGHLLMHLIGGATLALINISITTALITFIILVLLTALEFAVAIIQAYVFTLLVSLYLHDNT**PTKPTHTI**TQVPDLLQEPSQPY**RRA*PYDSTLTP--LSY*RQD*LPIS*QYISDDEM*SEKAPFKATTHQSYKKDFATE*SYLLSPKSYF-SQASSEPFTTQASLLLLN*ADVDHPQASTL*TH*RCHF*TPPFY*PLVSPLPEPTTV**KAIENKYSKPSSSQSP*VCTSHYCKLQNTMKPPLQSQMGFTAQLSL*PQAFMDYM*LLAPLS*LYASYAN*NSTSRQITTLASRPPPDTDTS*M*SDYSSTCPFIDEVHSSFSIKTVQLTSNQLASVHSGKEQLT***HY*QTPH*PLYWSSLPSDSHN*TPTQKKQAPM-NADLTL*DQPAYLSL*NSS**PSHSFFST*RSPSYFLSHGQP-KQQT*KPYSL*PSP*SHS*QSA*PTNELKRD*NEPN-MVFSLKQNK*FRLIKL*TNS*IPSVFSIYKYYHGLYNI-PR-RTVNISI-PLNILTPMSRRNDIITIY-HSN-SHH-P-KCTL-HPSQHNANYSTSFRSM*S-SPRTIATSNGIKHIRYRLRTKPKPSPMLKYIIPTIILIPLT*ISKNSII*TNTTAHSLLISFTSLLLLNQFNDNSLNFSPMFFSDPLSTPLLILTIWLLPLILIASQ-SHLLK-EPPTRKKL-FITILVT-LQTFLIITFSAIELILFYILFEATLIPTL-I-IITR*GNQT-ERLNAGLYFLFYTLIGSLPLLVALIYIQNITGSLNFLILQY*TQAVSNS*SNVFL*LACI-IAFMVKIPLYGLHL*LPKAHVEAPIAGSIVLAAILLKLGGYGILRITTILNPLTEIMAYPFIILSL*GIIMTSSICLRQTDLKSLIAYSSVSHIALVIVAILIQTP*SYIGATALIIAHGLTSS-ILFCLANSNYE-RIHSRT--IILARGLQTLLPLIAA**LLASLTNLALPPSINLVGELLVIMSSFS*SNITIILIGT-NIIITALYSLYILTTTQRGKYTHHINNITPSFTRENALI-A-L-HILPLLLLSLNPKIILGPLYCKHSLRKTLDCESTNRSSNPSYLPRKHARTANSCPHI*QYGFPRLLKDGSYPLVLGTKKLVQLQIKAINLFSS-TTLTILFVLTLPIIITNTNIYKSDKYPTYVKNTVSSAFLISLVPIIAFTNTGQEIIISN*H*ITIQTLKLTLSFKADYFSIVFAPVALFVTWSIMEFSIWYIHSDPHINQFFKYLLLFLITIIVLVTANNLFQLFIG*EGVGIMSFLLIG**HGRTDANTAAIQAILYNRIGDVGFIIAIA*FLSNLNT*DIQQIFIINPTHSNLPLIGLILAATGKSAQFGLHP*LPSAMEGPTPVSALLHSSTIVVAGVFLLIRFYPMIENNNLIQTITICLGAITTLFTAICALTQNDIKKIIAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFICSGSIIHNLNNEQDIRKIGGLFKTIPFTTTTLIVGSIALTGVPFLTGFYSKDLIIEAANTSYSNA*ALLITLVATSLTAVYSTRIIFFALLGHPRFPTSTLINENNPLLLNSLKRLIAGSIFAGFILSHNLPPITTPL-ITIPPYLKMTALAVTILGFTLAFEITLNTQNLKHKHPTNSFKFSTLLGYFPTIIHRLPPHLSLTASQKLASSLLDSA*LENILPKSIAHAQLKLSTLVSNQKGLIKIYFLSFLITIPLSLILFNPHA*LP*SPQHQQTRTNPLQQPIKHH-SY-IKLQHP*PPH*KP-QN-H-SY-HKQPNHPDH*T*TQAPPLL-PSRHIKPL-K-TLSSTLKETPPTQPY*RLKPRDTAQ*PLQSYNQIPPTSHPNKPKTPLNPKKTRQNSTPYHNQPHH*QSNPAHHKSAKALKKTQ*NQ-SQRQYSK*IQHTLSLFPYGL*P*PTA*KTIVVIQL*EH**QTSENLTP**KLSTMHSLTSQLHQTSHRDETSAPYLASA*SYKF*QAYSWPYTTHQIHSPHSHR*PTSAVM*TTGESSATYTQTAHPSSSSASLLT*DAAYTMAPTHS*KPETSELSYYSQP*LPRL*ATYCHEDKCHSEGQQSLPTYC-QLSPILEQT**NESEEAFP*TKPPLHDSLPSTLFFHS-LS-QH*PSSIYYSSMKQDPTTQQESPQTQTKSHSTPITQSRTS*VSYS**QHYS-H*PYLPQTS*GTQTTTPPQTPLAHHHTLNQNDISCSRTRFSDQSPTN*EAS*P*LSQS*SWP*SQYYTHPNNEA*YFDPSANACFEH*SPTY-*HSHELEDNP-SNTPSSSSDKSPQSYISS*S*Y*CP*QALSKTNS*NEESL*YMTLPRSCKPKKEATHLPETQGRSSSSTISTQS*NSR*TIP*FPCMYYLQDYKVLT*YYNLKCTY--IHGYVRRALLLYHIQ*YILCIIVHRTYYV*SCITLSSTCL*ACTSELLAPHGTC*FFIVHGTCTQIISSHQAYPAP*ITSLITRPRETSNP-LGRFPSSRSGPIACGGF*E*TLSGIWFLLQGHLI*NRSFFPLK*DISMD**LISPCSRITVISCIWYFYFLGDAWTQL*PPGLNLVN*LVAGL*LNLIYQPGSELQGLFSQWLQDIEKCDVHSFATGHAFIHCSLALARHAYYSVNGYRTYQYISPPCSFKSLIHLLKYASPKQRAIPLDSTNEFAKN*ISKPAQTTQGSPTS*T
    :: : ::  :|     |    |   | | |   |    :  ::    | : |:| :         ||   : :  ||    |: :    |:    |            :::   | :|:|: |  :    : :|:   :     | |   | :::   :    :|:|    |   ||  :|  ||     | |   : :  | |  | |   | | | |   :|       | | ::      ||  |     ||  |    ::   |  :     ::| | |:   |   |: :   | :|  :    |||      ||      : |:| | |    | ||: | |:|||||||||||||||||: : :|||:||:| |||||||||||||||||||||||||||:||: ||: |:|  : ::   | :|| :|       :  | | :: :  |   :| :  :: ||:|  | ||:|| |||||||:|:|::||||||||:   : | : | |::  :    | : | |  | :  : ::  | : : ::::  |     |: |:  |||   : :: :|  :|  |    :| |||:||| |:||:   ||  |    :| |||:||: | ||: ||  ||||||| |: : :  :  :  |    : ::  :       ::    :: |  ::| :   |: :  |: : :  | :|   | ||:  || |||||||| ||||| :|  :|    |:| |    |    ::   | :: |::   |||:|  :: ||:|:     | :::   ||||||||: ||:||||||||| :|    |      :   : ||::   |::      :       |:  |     :|   :   :  : |    |  :: :||     :  :| |:  : : :  :  :  : :|    ::|| |:: :    :: :::   | || :::| | |  || |   || :  |   |  |    :::     | : :        :  :: :|  : |:  ::  | : :|: :  |: :      |  |: :|::        | :   :| ::|   :::    |:   :| | : |  ::||  :   |  :  :|        :     :    | | :: |:::    | ::  ::  :  : : ::| | ||  :     |: ||      |  | :|:: ::| | : |    :| :: :::||:  || :  : :  : | : |: :: |   :    :: :|  |:  |  :::    : |    : : | :   : ::: |  : : : :  ::: : :  |     |  ::  : ::::   || |::|: :: :     |   |:    | |  :::  : :    |::   :    : :|:    ::|: | : :   :  :  : ||  |: :    |   : :| |  |   |: : |   :: |  :  : |  ::     : | |         ||    |  : : : : |:|    ::::    : :   |:| :   : : | : ::  |:|    : |  :  | :|   |: :: :    :|:     ||      | :   :|: : :  ::|: |    :  :::  |: :| ::: | :  :  :  :|  : ::| :    |  : | | |:|   ::  :   ::|    |  :: :   | | |    |   : |:||  | |  |  :  ::  ::| :   |   :  | |  :|:  :    ::|  |  |:| :  |:    ||   |  : : :      ||: || | :: |  || :::  ||     | |   |:  |||||| : | |  :: |      | :: :: |   | : :|    |: ||| |||||:|||||||||||| | |  | | ::   : :  :|   ::    |:  :| | :    :| | : : :|       :    | :|::: :::| |  | :  :  :| :: :     ::  :  |  : :: ::   : :    : : : :::  || ::    :     |  : | || : |  |     : : :::   |  : |  :   : :: || :  |  |: : : |        ::: :|: ||   |  ||    :|: |:  :: |:    :|     |  |      | | |||    :||  :  : |    |:   |::           |   ::: |:  :  ::  ::   : |  :: | : : :  | :| :   |  | |   ::  |   |:      |:  :   |  |: : |    ::     |:  | :  |::     ::: : :|       ::    :   | | : : |:  |  : :    :|:   :|: : |  :  | :: |::  |::  | :| :  |   |     |::  ::   ::      |:| :      : |  |: |    :  |   |   ::| :  :  |::: :||:::| | :   |  |::::  : |      |:|     |::::  :|  :::  :  :    :  ::| : |  :|: |    :::   :  | : | :   |:    || :  : |:   ::: ||   :|  ::  : |::|:::|  : |   :  |  |   :|    :  ||  :   :    :  :    : |:  : |   |:   |     :   : ::    ::  |     | :  |:: ::     | |  :   ||:    | ::     :|  |:| |::  | :| :| |    | :  |:|:| :|||     : :| |    ||::||||||||||:||||||| ::|||:|| :|:|:   ||:|||:: |||||:|:|||:::||| ||:||:|||:|||:||::|||||||||||||||||||||::|||||||:||||||:| | :|||||||||||||||:|||||||||:|||||||||||||||||||||||||::|||: :|:: ::|| | ||:||| ||:|||:||||||||||||||||||||||:    |||: : |  | |||| ||  ||| | | ||   |||   | |  | | :|  || |||  |||||:  || | | ||||| |:::| ||  |||||| |  |  ||  :   |||  |  |: | || |:|       || | ||||: |:|  |  | |   |      | |:| : |     : |||  |||||  |:  ||| || | || |||  | :   ||||:|||||   ||  ||    | |||   :: :  | |||| |: :   ||  |  |  |:|| || |   |  | || | |  |: |: :|  ||   || || ||||||   |::  |: ||::  | :| || | :     || || ||  : | :: |:||:||| ::   |:    : | :   :  || |:   | |   ::      :: |  |         | ::  :   |: :   |  | | : |   |:   :   :|| ||  :| :: |    |   :| :  |  :||  | : :  |   :| :   |    :     ::::  || :   | | | :|  | ::|    | :| |     :: : |::    |:| | : | | | | :: | :  : : || |:| : |  | |      |   :| || :  | :|  : : | |   : :: |:   : :     :  ::::   |:  | :: :  :|: |       |  | :  ::|           |:          |  :|  :  : | :  ::      ::: |  |: :  ||:  : :      :   : | : :       :|:   : ||: :  |: : | ||::    |:|  ::  |  :    :::     | |    :|: |   :     :       | :   : :  |:|  :|:   :  | :: |   ||::  : : : |  |:     : :::|  :  ::::|  |   : :| ::| ||||| :||||||||:||:||  : :  :     : : |||||||:|:: :: |||| |  |  ||: | :|  |:::|| |||: |:|:|:||| |   |    ||:||||||| | || :|:|||| |||||:| ||||||||||:|||: ||:|||:|||||||||:|||||::|||||||||||||||||||:|||||:||::|:|||||||||||||||||:|||:|:|||: : |:|| ||| ::| : |  ||:||:|||:|||||:|||||||||:|||||||||||||||||||:||||||:|: ||: ||||:|:|||||||||:|:|||||||||||:||||||||||||||||||||||||||||||||||||||:|||||||||||||||||||||:||||:|:|:| :||:||:|:|||||||||| |||:|| ||:||||| |||:|||||::||||||:::| |:||||| | |||||| ||| :||| |||:||||::::|: | ::|:  ||| |||:||||||:||:  |::::  |::|| | |  :| ||::||::|:| ||  |:|:| :||:|   ||| :|||::|||:|:: |:: | ::|:|||:||:|||||:: : |:|:|::    :| |  : :  ||   :  ::: |  | :| |   :  | |    |: : : ::  :   :: | :|   :| : | :    :| ::  |  |  |:   :: :  ::: |:  |  : ::  | ::  | ::|:|| |:  ::   :| :::  :  |  : :  :: | :: |:|  |   |:   |    |  : ::  :    :::|  |  |::   :     |  :|    :   | :   ::: :  :  : |   |:  ::|   : :: :  : |       ||  |    |  :  |   |||  :|    |   |   |  |   : :   | :  ::||   |  : | :     ::|:  :   : |:| : ||:: ::  || |  | :| | |::   : :  || ::|     ||  :|:: ::    :|   :| :||     :|   |  | |    :|  |||:   ||   : |   :   |  ||    |  :: : |::    | |     : | ||   |:::: :  | :  |  |::  |  ||:: :    | : ::   :    |      |:|   :: :::||  : : :|| | |   |:    |  |   || : |  |  ::|:  |  ||   :  :  : :|  |   : : |  : |  :  :  |    :|  |  ||    ::| : | :   |   | |||  |     :|  | : : |   : |  |  : : : ::|   | :|:  | :  | :  ::::     |   :|| |   :: :::  ::|:  |:    :    | |      | |: :| ::  |::| :   |:  :  |    :| ||    : : |  : : ::  | ::   ||:    | : : |  | |  : : :  |    : : |   
DHRSITLLTTHGSSPCIWYFRLGGVHAIALRDAGAGAPYVAVSVFDSCLILLFIAPTFNITGEHTY*SVLIN*CL*DIIITIECLHSRFPHRHHNKKFPPXXXXXXSGHST*THLCQTPKTKNPNTSLTRFQILSFGGMHF*QSPPN*HIIFPSHSHTTNLINTTLAHPTQHTHTHRC*PHTPNQPNPKDTPHSLCSLPPQSNTLKMFRRAHITP*TNRFGPSLSISS**DYTCKHPRSSEFTL*ITTIKRDKHQARSNAAQNA*PSHTPTGNSSD*PLAINESLTKLY*PQGWSISCQPPRSHD*PKSIEAGVKSVLDHPLPNKAKTHLSCKKLQLTQNRLRKWL*HI*THNS*DPNWD*IPHYA*P*TSTVKSTKLLARTLRATA*NSKDLAVLHTPLEEPVL*SINPDQPHH-LLL-SLYTAIFSKP**RLQSKRKYPRKDVRSRCSP*GGKKWA-TF---STPENYDSPYET*GSKVDLAVN*E*SA*LNRALKRVHTARHPPQV-YFKGHLTKT-PTHLYRGDKS*-HGK-CTGKCTWTNQSVA*HKAPNLHLGD-FNLT*PL*AKPSPKPTPPYYQTTLAKPFTQIKYRR*KLKP-GAIDIVPQGKDEKL*PSII*QGLTPIPSA**IN*K*LCKESQS*DPRNQTSYLRTAKRAHPSM*QNSGKIYR*RRQTYRAW**LVVQDR-ILVQL*ICPQNPLNPLVNLTVSPKRNSSLDTRKKPCRESKKFNT-HSRPKSSHQLRKRSSSTPTT*KIPNI*LNSSHP-IGPIYHPIEELMLV*VT*KHS-PPHKPASD*NTELTINSPISTINQQVIITLTVNPTQACS*GKVKKSKRNSANLTPPVYQKHHL*HHQY*RHRLPSD-TCLTAAVP*PCKGSIITCSLNRDLYEWLHEG-SAVSY-F*PVKLTCP*RGGHDTARREDPMEL*FINANST*QTHRS*TTKPALKISVG-ATSEQNPTSEQYMLRLHQS-KRTTILN-*S-NNLTNG-T--SYPR---DNSAI-L-F*S--PYQQ*GLRPRCWIRTSRWCS-RY*RFVCSTIKVLR-DLSSDRSNP-GRF-LST-S--NSSLYERTREIR-PT-SQSAFPRK*YHLNLVLYPHPPKNRVC*DGRAR*SHKT*NFTVRGSIPLLNNIP-MA---NLLLLIVPILIAMAFLMLTERKILGYIQ-L-RKGPNVVGPYGL-LQPFADAIKLFTKEPLKPATSTITLYITAPTLALTIALLL*TPLPIPNPLVNLNLGL-LFI-LATSSLAVYSIL*SG*ASNSNYALIGALRAVAQTIS-YEVTLAIILLSTLLISGSFNLSTLITTQEHL*LLLP-S*PL-AII*FISTL-AETNRTPF-D-LAE-GESELVSGFNIEYAAGPFALFFIAEYTNIIIINTLTTTIF-LGTTYDALSPELYTTYFVTKTL-LLTSLFL*IRTAYPRFRYDQLIHLL*KNFLPLTLALLI*YVSIPITISSIPPQT*EICLIKELL*WSK**ELKP-PYF*DYENR-THP*ESKI-LRATYHTPS*SKVS*ISYRAHTPKMLVIPFPY*LIPWPNPSSTLPSLQAHSSQR*ARTDFLPE*A*K*TC*LLFQF*PK-K*TLVPQKLPSSISSRKQPHP*SF**LSSSTIYSPDN-EP*PILPINTHH**S*WL*Q*N*E*PPFTSESQRLPKAPL*HPACFFSHDKN*PPSQSYTKSLPH*T*AF-SSLSQSYPS*QAVEVD*TKPSYAKS*HTPQLPT*DE**QFYRTTLT*PFLI*LFILS*LLPHSYYST*TPAPR-PYYYLAPETS*HD*HP*FHPPSSP*EACPR*PAFCPNGPL-SKNS-QKTIASSSPPS*PPSPSLTSTSTYA*STPPQSHYSPYLTT*K*NDSLNIQNPPHSSPHSSPLPRYSYLSPLLY**SYRNLG*IQ-TKSLQSPQ*VAIL-NFCNS*GLQNPTLHQLNANQPL*LS*ALTRPMGLKPTNT*LTAKHPNQLASIYFSRRRE-KRREKPRQV*SCF-FEFAIQYENHLGAGKKRPNPCL*IYS-PML-HSAILPHPH*CSPTVDYSLQTTKTLEHYTYYSAHELE-S*AQL*ASLFEPSWASQATF*VTTTSTTLSSQPMHL**SSS**YPS*SEALATD*FP**SVPPIWRFPA*TT*ASDSYLPLSYSCSH-LL*WKPEQEQVEQ-S-TL-P*QG--TTPTLEPP*T*PSSPYT*QVSPL--S*GPSISSQQLSI*NP---LP*PNTKRPSSSDPS*SQQSYFSYLSQS*LLASLYY*QTATSTP-PSSTPPE-EETPFYTNTYSDFSVTLKFIFLSYQASE*SPIL*LTTPEKKNHLDT*VWSEL*YQLAS*GLSCEHTIYLQ*E*T*TH-EHISPPLP*SSLSPPASKYLADSPHSTEAI*NDLL-QCSEP*DSSFFSP*VA*LALY*QTHH*TSYYTTRTTL*LTSTMSYQ*ELYLPS*EASFTDFPYSQATP*TKPTPKSISLSYSSA*I*LSSHNTFSAYPECPDVTRTTPMHTPHETSYHL*AHSFL*QQ*Y**FS*FEKPSLRSEKS***KNPP-*TWSDYMD-APHPTTHSKNPYT*NLDKKGRNRTPQSWFQANPMASMTFSKRY*KNHFITLSKLNYRLNPIYLNGTCSASRSTRRYFPYH--RR-AYHLS*SRPH-NHFPYLLPSPVCPFP-NTH-NKTN*Y*H-LRRSGNRNRLNYPARHHPSPH-RPPIPTHPLHNRRGQRSLPYHQINWPPMVLNL-RVHRLRRTNLQLLHTSPIIPRTRRPATP*R*QSSSTPD*SPHSYNNYITRRLALMSCPHIRLKNRCNSRTSKPNHFHRYTTGGILRSML*NLWSKPQFHAHRPRIN-SPKN-L*NRARIYPIAPPLPPLEPTVKLT*H*PFKLKIKRTNTSLQ*NAPTKYYRMAHHNYP-HTPYTIPHHPTKNIKHKLPPTSLTKAHKNKKL*QTLRTKMNENLFASFIAPTILGLPAAVLIILFPPLLIPTSKYLINNRLITTQQ*LIKLTSKQMITIHNTKGRT*SLILVSLIIFIATTNLLGLLPHSFTPTTQLSINLAMAIPL*AGAVIIGFRSKIKNALAHFLPQGTPTPLIPILVIIETISLLIQPIALAVRLTANITAGHLLMHLIGSTTLAISTINLPSTLIIFTILILLTILEIAVALIQAYVFTLLVSLYLHDNT**PTNHMPII**NPAHDP*QGPSQPS**PPA*PCDFTSTP*RSSY*AY-*-PTH*PYTNDGAM*HEKAHTKATTHHLSKKAFDTG*SYLL-PQKFFSSQDFSEPFTTPA*PLPPN*EGTGPQQASPR*IP*KSHS*THPYYSHQEYQSPELTIV**KTTETK*FKHCSLQFYWVSILPSYKPQSTSSLPSPFPTASTAQHFL*PQASTDFTSLLAQLSSLSASSAN*YFTLHPNITLASKPPPDTGIL*MWFDYFCMSPSIDE-GLTLLV*IVPLTSN*LVLTTFKKE**TSP*F**STPS*PYY**LLHFDYHNSTAT-*KNPPLTSAASTLYPPPASLSP*NSS**LLPSYYLI*KLPSFYPYH-EPYKQLTCH**LCHPSY*SSS*P*VWPMSDYKKD*TEPNWYIV*TKRMISTH*IMIIIFTKCPS-FT*ILY*HLPSHF*EY*YIAHTSCPPYYA*KE*YYRCSL*LLS*PSTPTPS*PILCLLPY*SLPPAKQRWA-*PY*SQSPTHMA*TTYIT-*TYS-NAKTN-RPN-NYITTTDMTFQKT---HN-LNQHNHPQPNY*HHPSTIF*PNQQQPI*L-F-PNLFLRPPNNP-PPNTNYLTPTPHNHGKPTPLIQ*TTITKKTLPLYTNLPTNLLNYNIHSHRT-NHIL-YLLRNHT-YPHLGYHHPMRQPARTPERRHI-LPILHPSRLPS-PTHRTNLHSQHPRLTKHSTTHPHCPRTIKLL-SQQ-LNMTSLHNSFYSKDTSLRTPLMTP*SPCRSPHRWVNSTCRSTLETRRLWHNTPHTHSQPPDKTHSLP-LPCTIPMRHNYNKLHLPTTN-RPKIAHCILFNQPHSPRSNSHSHP-NPLKLHRRSHSHNRPRTYILITILPSKLKLRTHSQSHHNPLSRTSNS-TPTNSFLMTFSKPR*PRLTPHY*PTGRTLCASNHVLLIKYHSPTYRTQHTSHSPILPLHIYHNTMGLTHPPH*QHKTLIHTRKHPHVHTPIPHSPPIPQPRHHYRVFLL*I*FNQNIR--L*I*Q-QRLTTP-YLPRKLTRTANSCPHV*QHGFLNF*RITAIHWS*APRILVQLQIKVITMHTTITTLT-LTSLIPPILTTLVNPNKKNSYPHYVKSIVASTFIISLFPTTIFMCLDQEVIISN*H*ATTQTTQLSLSFKLDYFSIIFIPVALFVTWSIIEFSL*YINSDPNINQFFKYLLIFLITILILVTANNLFQLFIG*EGVGIISFLLIS**YARADANTAAIQAILYNRIGDIGFILALA*FILHSNS*DPQQIALLNANPSLTPLLGLLLAAAGKSAQLGLHP*LPSAIEGPTPVSALLHSSTIVVAGIFLLIRFHPLAENSPLIQTLTLCLGAITTLFAAVCALTQNDIKKIVAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFMCSGSIIHNLNNEQDIRKIGGLLKTIPLTSTSLTIGSLALAGIPFLTGFYSKDHIIETANISYTNA*ALSITLIATSLTSAYSTRIILLTLTGQPRFPTLTNINENNPTLLNPIKRLAAGSLFAGFLITNNISP-ASPFQTTIPLYLKLTALAVTFLGLLTALDLNYLTNKLKIKSPLCTFYFSNILGFYPSITHRTIPYLGLLTSQNLPLLLLDLT*LEKLLPKTISQHQISTSIITSTQKGIIKLYFLSFFFPLILTLLLIT-*P-IPPSNFNYNIYTNKQCSTSNYY*STPIIIQSPRTNRILPNQP*PLSFINYSASYTIKVYHNHHPIILFHPQHQSYLHR*PH*NTHQDLNP*PPCLRIL-LNSHRCSISKDNHHSP*IN*KNY*THITSPKIQNNNTPDHTANNQ-Y*TPINRRRLRR-KPHKPHY*THTQQKQSIHH-YSRTDYNHDQ*-YE-KPSLYFNYKNTNDPNTQN*PP-NKIN*PLTHRPP-HPIQHLRMMKL-RL-T--PWRLPDPPN-HHRTIPSHALLTRR--LNRLFINRPHHSRRKLWLNHPLPSRQWRL--NILYLPLPTHRARP-I--LRIISL--LRNLKHRH-YPPAC-NY-SNSLHRLCPPVRPNIILR--GHSNYKLT-IRHPIHWDRPSSMNLRRLLSRQSHPHTILYLSLHL--ALHYCSP-SST---PPPIL-ARNG-IKQPPRNHLPFR*NHLPPLL-HNQ-RRPRLTSLP---SL-LNDINTILTR-PPRRP-R-QLY-P--SQ----P-LKHPS-PHQARMIF-P-IR--L-HNSPIRP*QTRR-RPCPITIHPHPSNNPHP-PYIQTTKHNISPTKPITLL-T---PS--RRP-PHSNLNRRTTSKLPFYH-HWTSS-IRTILHNN----P-NPNTNYLPN*KQNT-QMDLSL*YKLIHQSCKPEMKTF-FQGQIREKVFNSTISTQS*DSNLNYSL-FFHGE--ADLGTT-QVLTHPSTTAMYFVHYCQP-P*I--LYGTI-NT*PPVVHKNPIHIKTPSPC--L-Q-ASTAINPQL-SHINCNSKATP-H-PLGYQ--QTYP-P-LTVHST*SHLPYIAH-YSQI--PS-R-PH----G-*-----PP-SD-RGPLTTILREI-NIP-HK--S--A----TL-LAPGP*H-L-GV-A-K-VN--C-IR---HLV----PTS-G-P-*SLN--SPH-V--P---LK*DIT--M

showalignment(globalAlignment)
[score, globalAlignment] = nwalign(hippoPORF, humanPORF)

score =

 -102.6667


globalAlignment =

M-N------ENLFA-SF--ITP-TILG---L---P-LV-TL-IIM-F-PSI-LF---P-AP-TRLITNRLVSIQQ*L-IQ-LV-SK-QIMN-IHN-HKGQT*TLILISLILF-I---G-STNLLGLLPHSFTPTTQ-LSINLGIAIPL*AGTVIIG-FRNKTKISLAHFLPQGTPT-PLIPMLVIIETISLFIQ-PIALA---VRL-TA-NI--TAGHLLMHLIGGAT-LALIN-I--SI---TT--AL---I--T-FI---ILV--L-L---TAL-EF---A-VAII-QAYV-FT--LL---VSL-Y-LHD---NT
| |        |:| :|  :|   |||   |   | :|    ::: |  :| ||   |  | |  ||  :::    | |  |: :   | | : | : |    |   || :: |   | ::|    |  ::  ::| :| :: :|| | :  :| | |  :| |:  : |    |: ||  :: :| |::   : |: ||    :| :: ||  :|| : : :|:  | : :|| :  :|   ||  ||   :  | |:   :|:  | |   ||  :|     : :: : :: :|  ||   ||:   : :   :|
MANLLLLIVPILIAMAFLMLTERKILGYIQLRKGPNVVGPYGLLQPFADAIKLFTKEPLKPATSTITLYITAPTLALTIALLL*TPLPIPNPLVNLNLGLLFILATSSLAVYSIL*SG*ASNSNYALIGALRAVAQTISYEVTLAIILLSTLLISGSFNLSTLITTQEHL*LLLPS*PL-AII*FISTLAETNRTPFDLAEGESELVSGFNIEYAAGPFALFFIAEYTNIIIINTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFL*IRTAYPRFRYDQLIHLL*KNFLPLTLALLI*YVSIPITISSIPPQT

showalignment(globalAlignment)

>> [score, alignment] = nwalign(hippoPORF,humanPORF);
>> score

score =

 -102.6667

>> alignment

alignment =

  3×317 char array

    'M-N------ENLFA-SF--ITP-TILG---L---P-LV-TL-IIM-F-PSI-LF---P-AP-TRLITNRLVSIQQ*L-IQ-LV-SK-QIMN-IHN-HKGQT*TLILISLILF-I---G-STNLLGLLPHSFTPTTQ-LSINLGIAIPL*AGTVIIG-FRNKTKISLAHFLPQGTPT-PLIPMLVIIETISLFIQ-PIALA---VRL-TA-NI--TAGHLLMHLIGGAT-LALIN-I--SI---TT--AL---I--T-FI---ILV--L-L---TAL-EF---A-VAII-QAYV-FT--LL---VSL-Y-LHD---NT'
    '| |        |:| :|  :|   |||   |   | :|    ::: |  :| ||   |  | |  ||  :::    | |  |: :   | | : | : |    |   || :: |   | ::|    |  ::  ::| :| :: :|| | :  :| | |  :| |:  : |    |: ||  :: :| |::   : |: ||    :| :: ||  :|| : : :|:  | : :|| :  :|   ||  ||   :  | |:   :|:  | |   ||  :|     : :: : :: :|  ||   ||:   : :   :|'
    'MANLLLLIVPILIAMAFLMLTERKILGYIQLRKGPNVVGPYGLLQPFADAIKLFTKEPLKPATSTITLYITAPTLALTIALLL*TPLPIPNPLVNLNLGLLFILATSSLAVYSIL*SG*ASNSNYALIGALRAVAQTISYEVTLAIILLSTLLISGSFNLSTLITTQEHL*LLLPS*PL-AII*FISTLAETNRTPFDLAEGESELVSGFNIEYAAGPFALFFIAEYTNIIIINTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFL*IRTAYPRFRYDQLIHLL*KNFLPLTLALLI*YVSIPITISSIPPQT'

>> showalignment(alignment)
>> showalignment(alignment)
>> 
>> 
>> [score, localAlignment] = swalign(hippoProtein, humanProtein)

score =

   1.4023e+03


localAlignment =

  3×1967 char array

    '*HATTRHINMIYH-HPIHISDPIYYLSTENLKTHLPPKP*DYSSHNTKTAYPLRNEMNENLFASFITPTILGLPLVTLIIMFPSILFPAPTRLITNRLVSIQQ*LIQLVSKQIMNIHNHKGQT*TLILISLILFIGSTNLLGLLPHSFTPTTQLSINLGIAIPL*AGTVIIGFRNKTKISLAHFLPQGTPTPLIPMLVIIETISLFIQPIALAVRLTANITAGHLLMHLIGGATLALINISITTALITFIILVLLTALEFAVAIIQAYVFTLLVSLYLHDNT**PTKPTHTI**TQVPDLLQEPSQPY**RRA*PYDSTLTP--LSY*RQD*LPIS*QYISDDEM*SEKAPFKATTHQSYKKDFATE*SYLLSPKSYF-SQASSEPFTTQASLLLLN*ADVDHPQASTL*TH*RCHF*TPPFY*PLVSPLPEPTTV**KAIENKYSKPSSSQSP*VCTSHYCKLQNTMKPPLQSQMGFTAQLSL*PQAFMDYM*LLAPLS*LYASYAN*NSTSRQITTLASRPPPDTDTS*M*SDYSSTCPFIDEVHSSFSIKTVQLTSNQLASVHSGKEQLT***HY*QTPH*PLYWSSLPSDSHN*TPTQKKQAPM-NADLTL*DQPAYLSL*NSS**PSHSFFST*RSPSYFLSHGQP-KQQT*KPYSL*PSP*SHS*QSA*PTNELKRD*NEPN-MVFSLKQNK*FRLIKL*TNS*IPSVFSIYKYYHGLYNI-PR-RTVNISI-PLNILTPMSRRNDIITIY-HSN-SHH-P-KCTL-HPSQHNANYSTSFRSM*S-SPRTIATSNGIKHIRYRLRTKPKPSPMLKYIIPTIILIPLT*ISKNSII*TNTTAHSLLISFTSLLLLNQFNDNSLNFSPMFFSDPLSTPLLILTIWLLPLILIASQ-SHLLK-EPPTRKKL-FITILVT-LQTFLIITFSAIELILFYILFEATLIPTL-I-IITR*GNQT-ERLNAGLYFLFYTLIGSLPLLVALIYIQNITGSLNFLILQY*TQAVSNS*SNVFL*LACI-IAFMVKIPLYGLHL*LPKAHVEAPIAGSIVLAAILLKLGGYGILRITTILNPLTEIMAYPFIILSL*GIIMTSSICLRQTDLKSLIAYSSVSHIALVIVAILIQTP*SYIGATALIIAHGLTSS-ILFCLANSNYE-RIHSRT--IILARGLQTLLPLIAA**LLASLTNLALPPSINLVGELLVIMSSFS*SNITIILIGT-NIIITALYSLYILTTTQRGKYTHHINNITPSFTRENALI-A-L-HILPLLLLSLNPKIILGPLYCKHSLRKTLDCESTNRSSNPSYLPRKHARTANSCPHI*QYGFPRLLKDGSYPLVLGTKKLVQLQIKAINLFSS-TTLTILFVLTLPIIITNTNIYKSDKYPTYVKNTVSSAFLISLVPIIAFTNTGQEIIISN*H*ITIQTLKLTLSFKADYFSIVFAPVALFVTWSIMEFSIWYIHSDPHINQFFKYLLLFLITIIVLVTANNLFQLFIG*EGVGIMSFLLIG**HGRTDANTAAIQAILYNRIGDVGFIIAIA*FLSNLNT*DIQQIFIINPTHSNLPLIGLILAATGKSAQFGLHP*LPSAMEGPTPVSALLHSSTIVVAGVFLLIRFYPMIENNNLIQTITICLGAITTLFTAICALTQNDIKKIIAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFICSGSIIHNLNNEQDIRKIGGLFKTIPFTTTTLIVGSIALTGVPFLTGFYSKDLIIEAANTSYSNA*ALLITLVATSLTAVYSTRIIFFALLGHPRFPTSTLINENNPLLLNSLKRLIAGSIFAGFILSHNLPPITTPL-ITIPPYLKMTALAVTILGFTLAFEITLNTQNLKHKHPTNSFKFSTLLGYFPTIIHRLPPHLSLTASQKLASSLLDSA*LENILPKSIAHAQLKLSTLVSNQKGLIKIYFLSFLITIPLSLILFN-P'
    '|:| |::  | :| :| |    | :  |:|:| :|||     : :| |    ||::||||||||||:||||||| ::|||:|| :|:|:   ||:|||:: |||||:|:|||:::||| ||:||:|||:|||:||::|||||||||||||||||||||::|||||||:||||||:| | :|||||||||||||||:|||||||||:|||||||||||||||||||||||||::|||: :|:: ::|| | ||:||| ||:|||:||||||||||||||||||||||:    |||: : |  | |||| ||  ||| | | ||   |||   | |  | | :|  || |||  |||||:  || | | ||||| |:::| ||  |||||| |  |  ||  :   |||  |  |: | || |:|       || | ||||: |:|  |  | |   |      | |:| : |     : |||  |||||  |:  ||| || | || |||  | :   ||||:|||||   ||  ||    | |||   :: :  | |||| |: :   ||  |  |  |:|| || |   |  | || | |  |: |: :|  ||   || || ||||||   |::  |: ||::  | :| || | :     || || ||  : | :: |:||:||| ::   |:    : | :   :  || |:   | |   ::      :: |  |         | ::  :   |: :   |  | | : |   |:   :   :|| ||  :| :: |    |   :| :  |  :||  | : :  |   :| :   |    :     ::::  || :   | | | :|  | ::|    | :| |     :: : |::    |:| | : | | | | :: | :  : : || |:| : |  | |      |   :| || :  | :|  : : | |   : :: |:   : :     :  ::::   |:  | :: :  :|: |       |  | :  ::|           |:          |  :|  :  : | :  ::      ::: |  |: :  ||:  : :      :   : | : :       :|:   : ||: :  |: : | ||::    |:|  ::  |  :    :::     | |    :|: |   :     :       | :   : :  |:|  :|:   :  | :: |   ||::  : : : |  |:     : :::|  :  ::::|  |   : :| ::| ||||| :||||||||:||:||  : :  :     : : |||||||:|:: :: |||| |  |  ||: | :|  |:::|| |||: |:|:|:||| |   |    ||:||||||| | || :|:|||| |||||:| ||||||||||:|||: ||:|||:|||||||||:|||||::|||||||||||||||||||:|||||:||::|:|||||||||||||||||:|||:|:|||: : |:|| ||| ::| : |  ||:||:|||:|||||:|||||||||:|||||||||||||||||||:||||||:|: ||: ||||:|:|||||||||:|:|||||||||||:||||||||||||||||||||||||||||||||||||||:|||||||||||||||||||||:||||:|:|:| :||:||:|:|||||||||| |||:|| ||:||||| |||:|||||::||||||:::| |:||||| | |||||| ||| :||| |||:||||::::|: | ::|:  ||| |||:||||||:||:  |::::  |::|| | |  :| ||::||::|:| ||  |:|:| :||:|   ||| :|||::|||:|:: |:: | ::|:|||:||:|||||:: : |:|:|:: |'
    '*NAPTKYYRMAHHNYP-HTPYTIPHHPTKNIKHKLPPTSLTKAHKNKKL*QTLRTKMNENLFASFIAPTILGLPAAVLIILFPPLLIPTSKYLINNRLITTQQ*LIKLTSKQMITIHNTKGRT*SLILVSLIIFIATTNLLGLLPHSFTPTTQLSINLAMAIPL*AGAVIIGFRSKIKNALAHFLPQGTPTPLIPILVIIETISLLIQPIALAVRLTANITAGHLLMHLIGSTTLAISTINLPSTLIIFTILILLTILEIAVALIQAYVFTLLVSLYLHDNT**PTNHMPII**NPAHDP*QGPSQPS**PPA*PCDFTSTP*RSSY*AY-*-PTH*PYTNDGAM*HEKAHTKATTHHLSKKAFDTG*SYLL-PQKFFSSQDFSEPFTTPA*PLPPN*EGTGPQQASPR*IP*KSHS*THPYYSHQEYQSPELTIV**KTTETK*FKHCSLQFYWVSILPSYKPQSTSSLPSPFPTASTAQHFL*PQASTDFTSLLAQLSSLSASSAN*YFTLHPNITLASKPPPDTGIL*MWFDYFCMSPSIDE-GLTLLV*IVPLTSN*LVLTTFKKE**TSP*F**STPS*PYY**LLHFDYHNSTAT-*KNPPLTSAASTLYPPPASLSP*NSS**LLPSYYLI*KLPSFYPYH-EPYKQLTCH**LCHPSY*SSS*P*VWPMSDYKKD*TEPNWYIV*TKRMISTH*IMIIIFTKCPS-FT*ILY*HLPSHF*EY*YIAHTSCPPYYA*KE*YYRCSL*LLS*PSTPTPS*PILCLLPY*SLPPAKQRWA-*PY*SQSPTHMA*TTYIT-*TYS-NAKTN-RPN-NYITTTDMTFQKT---HN-LNQHNHPQPNY*HHPSTIF*PNQQQPI*L-F-PNLFLRPPNNP-PPNTNYLTPTPHNHGKPTPLIQ*TTITKKTLPLYTNLPTNLLNYNIHSHRT-NHIL-YLLRNHT-YPHLGYHHPMRQPARTPERRHI-LPILHPSRLPS-PTHRTNLHSQHPRLTKHSTTHPHCPRTIKLL-SQQ-LNMTSLHNSFYSKDTSLRTPLMTP*SPCRSPHRWVNSTCRSTLETRRLWHNTPHTHSQPPDKTHSLP-LPCTIPMRHNYNKLHLPTTN-RPKIAHCILFNQPHSPRSNSHSHP-NPLKLHRRSHSHNRPRTYILITILPSKLKLRTHSQSHHNPLSRTSNS-TPTNSFLMTFSKPR*PRLTPHY*PTGRTLCASNHVLLIKYHSPTYRTQHTSHSPILPLHIYHNTMGLTHPPH*QHKTLIHTRKHPHVHTPIPHSPPIPQPRHHYRVFLL*I*FNQNIR--L*I*Q-QRLTTP-YLPRKLTRTANSCPHV*QHGFLNF*RITAIHWS*APRILVQLQIKVITMHTTITTLT-LTSLIPPILTTLVNPNKKNSYPHYVKSIVASTFIISLFPTTIFMCLDQEVIISN*H*ATTQTTQLSLSFKLDYFSIIFIPVALFVTWSIIEFSL*YINSDPNINQFFKYLLIFLITILILVTANNLFQLFIG*EGVGIISFLLIS**YARADANTAAIQAILYNRIGDIGFILALA*FILHSNS*DPQQIALLNANPSLTPLLGLLLAAAGKSAQLGLHP*LPSAIEGPTPVSALLHSSTIVVAGIFLLIRFHPLAENSPLIQTLTLCLGAITTLFAAVCALTQNDIKKIVAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFMCSGSIIHNLNNEQDIRKIGGLLKTIPLTSTSLTIGSLALAGIPFLTGFYSKDHIIETANISYTNA*ALSITLIATSLTSAYSTRIILLTLTGQPRFPTLTNINENNPTLLNPIKRLAAGSLFAGFLITNNISP-ASPFQTTIPLYLKLTALAVTFLGLLTALDLNYLTNKLKIKSPLCTFYFSNILGFYPSITHRTIPYLGLLTSQNLPLLLLDLT*LEKLLPKTISQHQISTSIITSTQKGIIKLYFLSFFFPLILTLLLIT*P'

>> showalignment(localAlignment)
>> hippo_random = hippo(randperm(length(hippo)))

hippo_random =

    'CGGTAACCGAAATGCAAACCTAGACCCAACCTCACAAATTCCGAGAATCGATACTGTTACCTATTCATCCAGTTTTCTCCGCACTTTCGACTGTATAATCACACCAGCGGTACCCTACCGGACGCTCCTCGACAAATTATCTACTACAAGACTATTAATCCCACCGATTCATCCTTCTAAAAAACTGTGCTCTCTCGCCCAATCTAAAAACCTTACCTAACTTTTGAATAATATTCGACAGGTCGATAAAAAAAATGCTGATTGACTATTCTATACAATAGGACATCTCCCAGCCGCGAGTACTATCACTGATATCACCAAGTATCTCCTACTACTAAAACAATACTATACATCGGCTCACTAGAACAGAATAATTTAAAAAAATACTCAAGATGTACCGGCTCATCCGACTTGTATCCCAACCTTTTTACCATTTACTAATGTCAATGGTGGACGAACCTCATTGACGAATGCATCTAGCGGACTCATCTCCAATAGCGCCACAAACTCCATATGGTCCCTTCTCACACATTACACCACAGTAAACGATTCCCATTAGTAAACAGCCGATAGTGCTAATTCACCTTCTACTTCCCCTACATAATTATTATCTCCCAACCGGTACCTAATTGCTTCTGGCACATGAGTTCTCGATTACATCACATAATTGTCCCTATATGCGTGTCAAGATCATGCCGGAGATTAGTATGGGGATCGATTAACATCAATCAAACCTCAGGAATCAAAACTCCCTACACTAACTTAGAGACCTTAGCTTTTTAGTGACACAATGCAAAACAACTAATATACCGAATATCCGACCCCACTACCGCTACAATCTAGACATTAATGCACATACTAGTATCCTACAAATGAACAGCAAACCCATTGCTACAACGCTCACAAGAACTGTTTCTATAATGTCCGCCTTCCATCTTTCGTGAGGGAACACCATATCTGTAGAAAGATATAAACGCAAAATAGTCGTTTTACCCAGAGATACTCACGACTACGGAGGCCAACGCAACCCACGCCTTCAGCGATTGATGCCACTATACGGGTTGCATCGGTCACTTAACTTCACACTTCTAAAGTTCCGCCGCAACGCCTCTTCTCCATTGACTAAAGGTCCCCCACCGCAACTCAATAACGTACACCTATCCTATATCCGTAGAACAAACAGATCGTCCCTCTTTCATTACCAACCCCAATTGAACTAACAACATATAGCCCTAGGCCCACTATATCGCTGCAATACGTCCCTAAGCAGGCTAGTTCCAATAATTCAATCTAACTTAATGGCGGTAACACAAAGCACAATTCGAAAGTCAACAGAGCCCCGCAAACACTAATAATAGAACAGTATTTAATCACATATCAAGACCCCCCTAGCTCACTGCCGCACTAAAATAACGGATGATAATAATTGCCCCAGAGGCGCTTACTCCTAAACTAAGCCGTAAGAAAAAGTACTTAACTCCTGTAAGTCAAATCACGGCCGCAAATACATTCCCTCCCCCCTAACCGGTCTTGTCGCCTTTATAACGTAACTGCAAAAATCATCGCACCGGCTCCCTGTCTGGTCTCTTGAGCACGTTACCTCACCAGGACCTGATTTCAATTCACTGACTGCCTATCCGCCTCTTCTCCAAGTAGTATCTCCTAAAAAAATAGAGACGGTCCCCACCGCACTACCGTCCACTGACCTCACTGCCTAAATCCTAAGCACCAAAAAACTCGTTTACCAAGCTATGCTAGCTACCCTTTGTTGTACGAAAACCTATCTGTCGATGTCTAAACGACTGATTAAAAACAACATAAACACTCATCACTTTCGTAGATACAAAACGGATCTACCTAATCAAAGATATACGGCATATGCGGGTATGCGCTTGTCATGTAACAACTGCCAGCGCTTCAAGCGTCTGCGCACAGAGGGCTTCTCATAATTACGAAGCAAGTTACCTATGACCCTACTAGGTAATCGAAGGGACCACATTCTTCCCCACCGCTAAAGAAAACATGTAATAACCGAATCATCAAACACAAACAGACCTACCAATCCGTAACCCTCAAGACTACATAACCAACCGCTTAAAACAAATTAGTATATAAGCTAGAAAGGTGACCTGCATAAACTCGTCATCTCAGCAATCCCGAATAAAAACCGTCTACCATTGTACAAGCTAGCAGTCAGATTCGAGACAAGCATCCTCAAACTGTATCATAATTACAGATTAGCGCACAAGAGCACAAACGAAACTCGTCTCCTAAGCGGTGCCTTAAAGAATCAGCTAAATGACCAGCTACCAAACCCACACCTCCTCCATGCAAATAAGTCATAAACTATAGCCTGTATCGTGATGAACGTACTTACCGGTACAAACGATTTCGTACATACGTAACAGTAAATTCCATTCCCTTAATACCTAAACCATTCTCAATACAACCTCTTGCGGCAAACCTAAAGTGTCCACCTCCATTCGCCTTGTCTCACCAACCTTAATGGTTCCAATCAACAATCCTAACCTCCACTAAAGTAAATCTCAATCTTCCAGATCTACGCGATATGATAGGCAACGAATCTTCACACATTGGACATTCGCGCACTCTTATATTGCCTACCGACGTATTTAAAGTCCATGGATATCCACTTAAAACACCAATCCCTGCTCCCTAATCAGGCGCTAGAAAATAACCTCATAGACCTAGACACTCCTATCACAAATACGCATGAACGTTATTACCTGACCGCATTCCTGACCGGAAATCCATGATTTATACGTTCAAATCGTTCTCACTTGCCCAGTTACACTCCGTACTAGCTACGAGCCAAGAACTAAAACTCGGCAACCGACTCCTATAACTCCCTCTTCTTCTAGAGAACTCAGGTACGCATGGCTAAAAAACAATGCAATATCATTAAAAATACTCTTTAAGCCTGACAAAGAGAGACACAATCCAGGTGACTCGTTACTCCATGTGAGAACTATTAACACTTACCAAGGCTTAAAGTACTATCCACCTATTCCTCAGGAGCCGACCTAAGGTCCACCTGTATTAAATCACCCATTTGGGAACCCACCGGGAACTCGCAGTTTCAACATACCGCAAACAAAGATTCATAAAACAACCAATAACCTCCAGGAATCAATCCAATCACTGTACAAACAGTCCACACGAGGCGAGCAACTATTGATGGTATGACACTACATTTCGAGCATCGGGGCAACAGGCGAAAGTCTCATCATAAGAACACTTACTTTGCCACCACACATTGGAAGTAATTCTAAGACTCCTAGTAAAGGTAGGATCTGTTGTCAACCATTGACGTGACGTGAGTCAAAGTTCAACAACACCCACCCGTAACACAATCTCCTGTAAGCCTACCTGTCAGCTGGCAATCGATTCATCCCTAAATAGAAACGCTAGGAACTTCTCATCAGCTATAACACTTTCGAATGCACACCGGGTGACCTATATTGGTTAAGTACTCATTTCGCCACAACGCCTATAATCCGTAGTAATCTACGTAACTAGTCATCCATATACCGATACTGTTGTAGACACCCGACTCTCCCTCTCTCAGCCCTCCAGTCATGAAACTAAACTTTACGAGAAGCCGACTCAATTGTTCTTCTCCTACTCATGTATAATATTTCTAGATCACTCAGTCCCTTACAAGGTCCCCTTAAGAACAGAGCCATTCACACCTCAGCGCAAATCTTATGAGCCGATGGCATTTCACGTGCAGAGAACCCGTTACTGCAAATCCCGCCCAACTCCGCGGTAAGCATAAACCTCCAATCGTCCTCTCACAGCCCCTCCAAATCACTCCAAGCACGTTAAGCCAGTAACGACTTGTAACACTTCATCTCTTACATCGAACCATAAGAAATTCACATACCCTTGTCAATTACCGGAAAGACACGCACTTGATAGTTATATAGACTCTTAGACGAATACATTCATCACACATTTTTATTCCCTCGCAGATCAACCACTCACCCACACTCCCGATTAAAAATTCCATACCATCCAGAGACCGATCAACCTGACGTGATCCAATCAACCACACCCGTCCCTGACATCCCTATGTAGCACCGACCGGTCAACAATAGGATGTTGCTTATGGTTAAGGACACACTTGAATATCTGCCTCACAGAGGTATACTACTTGTCCATCGGATAATAGATGTATACGAATAAGAACTATCACCACTGCGACCAACACACAACACATCGATAGAAAATAAACACCAGACATACAAGCGGCGACAAATTAAGAATAGCCGACCCCGATAAAAACGGCGAAACCAAAATGCTTCTTCTCTGATTTATTTCGGCATCACATGACCCCTAAGCTCTAAACTCGCTATTAAAGTACCTATGAACTCCGCGAAAATGCTCCCGCCCCGCCCAAACAACCTCTTTCCTGCACCGGAGCCCCCCAATCTAGAAGTTACCGCTGGGTAGCGACTTTAGGGAACCCACCCACTTACCACACACATGCCATTCGTTGCCCACCGTCGTGTTTCCCAACCACCAATCTTCGTAAAGTGGAATATCGAAAGCGGCACTCTATTATTCCTCTGTTTCCGCATAAGCAGGAGGAAACCACGATAAACATAATAAATTTGTTGCACACCACCCACACAATACCCTTCTCACCCGACATTATATGCCAGTAACACCGCCCTAACAGACCCGGGCTAGAAAATACTTTCGTAGAAGCAGACGGCAGGACCTACACGGCCATCGCACCCTGAAACCCCCCAGCTGCTTTAGTACCTATATTCTTTGCGTCAACAACCAAACTATACGTACCTTCATTGCAGAATTTCCGACAGTGCAACCTGAGCAGTAAGAATCTTGTAACTCTACTGAAAACCCAAACGGCCATAAGCCGCTTCCCACAAGATTTTTTACTTCCAAAATTCGTCATCCCCTTACAAATTCCGCCGCACCTAACACACATAATACATCAAATCGAACCGTTAGCGACAAAAGTTGATCGACAGTAGAAATAAAAAGTACTTGCTCAACATAGATTCCACCCAGATAGCATCACAATCCAGTTCTGAAAAACCTTCACACCTACAGCTAAAAACTGGTGGTAGCTGTTGGATAACACCAAAAAATTTGACAGGTTTACCCGCGTAATACAGACAAGACATTTGGTAAACAGACTGACTAATGCTTTTCCAAACATAACAAACTAAGCTTGATGTAACATCTTTCTCAACCTATCCTTTCTAGCTTCTTCCAAGTAACAAATTTTGACCGATGCTATCTCCTCTTAACCAGGACCGAAACTCTTAGGAATATAAGCTGAGGAAAATAACAATCTTACAACGCATCTGCACCAAGGACCGAGCCCAGCGAGAAAATCGTCTTCAACACGCACCCCAATATCATAGACTTTATAGCAAGAGACCCTGCGCCCTCCCTTATTTATAACTTTACGATATCATAAAATCAAGTGTTATAAAACAGATACCATTTTTACCTCAATGTTTCAGACTCACTACACCACTAACGTTGCAATAGAAGATCTTGAAAAAACCTTGCTTCAACTGGATTATTATTCAAAATCAACATAAGCCCATTGGAATCGCCTACGAAGACATATCTCTATCCACCTCCGCCGAATGGAAGCCAACCTAAGTATCGCCACCTCATACCCCTAAAAATGTTAAACATTCCACCGAGGTTGCCCGCTCCATTCTAAGACTTTTAAGTTTAGGCCCATGGCAAGGTATGCTCGCCGCGCCACAACGAACTTGGACTATAACGGTCACCGACGTATATCAAGTCAGACCCAACGACTTCGCGCCCCGTATTGACAACGGAGGCAAATTCCACGAAGCTCAGCTGAGCACACAACAGTAAACGCATACGCTTTATTTCACCTACAAATTCCTGGTCACTCAGGATATCTAACCACATTAATGCCAAAAATACAACCATACTCCGACACGTGAAACTATCCACAACCCCCTCCTTTCACTAACCGCCCTTTAATACCGCTTTAAAGGCTACAATTAGACAATCTACAGCATCATTCAAATGAGCATTCTCCCTGTAATACTTTCGAGCTCAACTACCCCGTCGATCCGGAGTGGCAAATCCTTCCTTGCACTAATACCATTCAAATCAGATAACACAACTCCTATAGAACAACACTAAGACCTTTAACTTAAGCACTCATAACCTTAACCTTCACTGCTCCAAACTTCAAATCTCTTCCAGGACAGGAACTATTGTAAACGAGAGAGTTTAAAAATCCTATTATCTTAAAGCACGACATTATCCTGCGCACTGGAAGTATCCAACTTGGATCCCAAATTAAGATTCCATTCATGAAGAATACTATAAACCAACCACGAAGAACAATATAGGCATCTAAAACTCCCGGTCAACATAGCCGAAACGTCAAGGCCTAACAGACCTTCTCCGGAATTCTCTACTATACCAAAATCCCCTGATCTAAAAAATTAAACTGAGCTTAACAAAATGGCAATGTCTTAGACTTACCAACATCATACAACTTGACCAGAATACTCTTTTCCACAAAACCGATGATAATCATCGATCAGTCCACCCGTTAGACAAACTGAAACTCAATAAGAAGCAGACCGCACCCTACTAATGGCTAAATTATATACTAATACTTGATAACTCATCTTTTAACATCGAATCGACCCGTTGCAGCCCCAAGAATATCACTATCGCCAGTAGCTTCGCCAATTTCTTAACGTTGTAAAGTCCGAAATTAGTATTTAATGGCATCAAGGTTCGAAACACACTCGCTATCCGCGACAAAGTAAACTACTCACTAAATCCTCCAATTATGTCTTTCGACAAACCAGAATTTACGGTGAATCACTCAAGCCCTCCTCACCAACGCTAGATAAAATCGGCATCTAGACGGAGTCTCCGTTTAATAAAATAACGGAATTCGACGTTCCCCACACATACCAAGGGCAACGCTTTACTCCCAAACATTGAAGTAAAAATGGGGCTTCTTAGCATCCATAGCCAGGGCGAGTCCCCTAATACTCCTAGTGAGAAAGCGAAAGGTAAAGCGCGTCAATCTAGCGAGAATACCCAAACTGCTTAGTCACGATCTTTGTCCTACATAATTGACTCCTATGGCATCTTAGGAAAGCATCCGACAGAGGATCTATCCTGCTGATGCATTGCGTTCGTATTAACCAACAATTTACCTGTGAGCTCTAATCTTATCATATGCCCCAACAAAACTTAATCGGTTGCTACATCCAGCAATAATTTTCACCTCACATTCAACTGTCCACGTATGCATTAGCTCTATCCCATTACCCTAAAACCTATCGCCAATCTTAAGCGTCAGACTCACGTAACCTCAAAACTGACTACCTCATTCCAAATCCAAGACTACGATCCTCCCACACATGTTAGCCACGACGATAATAAGCGAACCAACCGTGGAAAATCACATTCGACTTTGTGCTCACCTCTACTGACCATCGACCACGTACACTTCTGGTCGTCAACTATGCGGGCCCTCCGCCAGACTTACAAGCAACTTCTTCGCAGAAACTGACATTAAGAAAAATTAAATACAGGCGCCATACGCCCCAAGTTATCCATTCTTAGATTCACTTTATTCCATCAAACTATGTATTCGCAACTCACAAGGCTCCCAATGCTGACTGACTACTACCCAGACCTTCATATAAACTACGTGCAACTGTAAAATGCTTCACGACCATTCGCTTTCCCCACAGTCCAGAACCTGACCCTAAGAAGGTAATATAACGCCAAACCATCGCTTCAACAAAAACATAATAGCAATTCCAAAATGAACAGGATAGATAAGTTAAAGGTTTACAGTCTTCCCCATAACTATGTTCCTGATCCTTAAATATAATGGGTTAGTACAACGCAAACGCTCAAATCGGATCATTCTTGGATAATAACAGTCCCATATACTTTTCAGAAATACACACGTGAGGGAAAACTAGAATGCCGTAATATAGACACTTTCTGTCGAGAAGAGGTGATTTCGTTCATCCATCCTTAAAGTTACCGCAACCGCTTCTCTAGCCGATTCACAGACCCTTCCCCCCAAGCTTACTACTCCCTGATCACTACCATTAGCTAGTGTACTGTCAGGCCGGATTCATATAGTCACCACCACGCTATGGTCACCCCTAACACTAACCTCCTTAAAATACCAGCAATAGGCCATACTCCACGGATAGATTACCTAAAAATCACACCTAACGCGAACAAACTAGCCTATTATTCTCACAATTGATTATTACACCAGACATCCTAAATTTCTTCTTGCTAATCCCTATCTCTCGACGCTAGCAAAGAGACTACTCCTAACGTTGTTTCTAAACACCTCGCCAATCTCTTAAACCGAAGTTTGACAATACCGTTTCTCAGAACCTCTTCTAAATATGGCAAGCTAACTTTCATCTTAGAGTTATCCAGGTAGTGAAAAGTATACGTTACGTGCTATATCAGTATCCGCAACTATATACTCAAGGATAATCTTCCAGCACACACTTATATACACCTTGATTTATGTATCATAAGAAATACGCCACAGTCTTGAACTGAATCGTGTGAAAATTCTACTCCATGAATCCAGTGACCTTACGAATATCGCGAAATCATCACTGAGCGTTCTTAGACCCTTCTGGCACCTTTATGTTTTTAAAATACTAAACATGCAAACAACTTCTTATTGTAAAGCCTAGCACTTCAGCCACGGATCCTCCGAAAACACCAACTAGATCTGAGTGACAACATGAACTCCCCTTAATAGCGATTGAGCACTCTCATTTACCGACGCAACTTTTGGATCAATTTATCCCCCCCTGATTCAATGACCTCCGCATCTAGCCACCGCTTCTAACATGTGTTAGTGTACAATAAAAGACCCGTACAACCCAAAAAATCTGCTCACTCCGCTATCAGCCCAAGGACATCAGCGTCATGGAACCCCATATCCAAACCTCAATATAACACACCCCTCGTCTCCGTAATAATTATAGTTACTCGGATTGTAAACATCCCTCGCGTATTCTTCTCAGTACAGTAATAGCTTATCTTCCCGATGTGCACCCTCCTGCGCGGTTTACCTAAACCATACTCATCTTCCACCAATCTGGTACATCAATTCCATAACACCCGATGAAGTTTTAACCGCCTGAGTGTTATTCAGCTTCAGCATGTCAAACTTACCAAATGCCGATTGTGGAGCTCACACATTATTGAACACAACTCTTACTGACTTAATCATTCTACCCTACCGCCATAAAGCCTGTTACGAATAGCGCATTTCTAACATAGAATCACCAACCATCTAACGATGATTCTCACCCCATACCCCATCATCTAAGTGCCACCCTCTAGCTAGAACAGTTCTTCCACCGCAACGTTCCGGGTCGTACACTTAACGTTACAAACGAGCCTAAAATTCAGGGTATACTCTGCCCTTTGAATAACTAGTCATATTCAGCTCTATAAACCAAAACATTATTTTAACTTTCCCTAAATTTACGTCAAATAAACATTCCTTATCGCTTACCTCTTTAGTTACTGCTCGACGTACGCATGCGACTACACTACGGGGTACCCAGAAAAAAGCTCAACCACGCAATAACAAACTGTGCATGAACCGTTAAACTTACTAGAAATAGGTCTCATCCCTGCCACTAGGATCCAGACAATCAGAATAGGACACTTCCCCTAAATCAATTGTAATCTCCAGAGGACTAGCCAAATCATGCCTCGTAAAACTAAATCATAGGCACCTCCCTCGCACCGATACTCTCCCTAAAAAACACTATCCAGATCCCCAACTACATAAAACTAAAATCAGACGCCTAAAACTAATATCTTCTCAACTATAAAACCGATGAACACACACCAATCACTACTTTCCTTGACATAGCCCTCACCAGGAATTTTTACGCAAATGAATATCTCTGTACCACTGACTTAGAACGATTCGACACAAGAACACAAACAAATTCATTTAAGTTTTTTCGCTTTACTTTCGTATACTGTGAAGACTCCAGCCCTGTAGACCTTTCCCAACCATTAAATCTAGTCTTACTAGACAGTATAAGACTATGGCATACGGTCCTATAAATACTAATTTCTGAACTCTATCAATTATACAACTCATATTATGAATGTCCTGACCAGTACCGAGCCAGGCGAGCGATGATCAGAGCGTTAAATATCTCAAAGCCATCCCAATTAACGCCCAAGAACCCTCACTAGAAATAGTGCCCCGAGTTCTCCCCACTATAACATCACCACTATACACCTAGGGATCATACATATAGACAGCGCTCAATGGTACTACTAAAAGTCGAAACGTTTTGAGAGAACAAAACCACATGTAATATCCCGTGGATTCTAAATTTCAGGCCAACAGTGTACACAACACGGAAACTACCACTAATATCTTAATAAAAACCCCAGAATTCAGGTCTAAGAAAACATAGATACCTCATCCTCTCCCTCCGTGCTTAAGGCAATTAGCATCCTCTGACATATAGAATTTAATATTACTAAGAAAACGACAACCCAGGCGATTCTATTGTACACTCGACTATACAAAGTGTTTAGCTGACTCTCGCGCACCCGTATCATTCTCTCATTCCCCTTTCTCAAAAACACGTGCCCTTAACCAGCCATACGACAAATTTTCCAACATTTAATCACAAAATGATCTAAAGACAACCATGCATTCGAAGGATAATCCCGCGCCCTCATTAACCTTAACGCTCTCTTACTGATCCCACGACAACGACACTCTTACTAAACTCAACTACAGCCTCTACGCTCCAAGTTAATCCTTACAAGTTTTACAAAGTACCGAAGGTACAATGACCCAGCCTTTGGTATTCCACTCAGACTCCTCCAAATACATATTTGGAACAAGCATCAGAAACCATTCCCTAAATCTTTCTACAGGGAAAGAACCCTAACAGGGAATCAAACACCTCGAACTTTGTCCTCTCTACCGACGGCACATCGAGCATCCCCGATACCATACAGCACTTGCCGCTCATACCAAAGAATCGCCGAACAATCGGGTCCCTAATTACTATTTAGATCCTCCAACTCTCTCTAACTCCCAATAGTAGCCTTAGTGTAACGACTCAAATAACCATTATATCAATAAAACAGTGCCATAGAAAACTATTATACCGCGTTTCAAAACCTAAAACCCTGGACATCCGCGCCCACAAGAAATTCCAGTCCATCCACCTAACAAGTCTACAGATAATAGCTACACATAACAGGCACGACATCAAAAACAACCATTATCCGTGCATCTCAGATATTCCAGCTATATCAAAAAATGAATTCCCGTTCGCACTAACTAAAACCACCTGTTTAATTATGTACGGGCTCCATTCACTTAACAACAGAGACAGAACTCTCCAGAAGGTCCTTGCTGGAATAACATTTCAAGTCCCTAACCATCCCAACGAAGATTCATCAGCACTCTGGCTGAGAATCTTCGAGACAGCATTTAATCGGCCGTCATTGCACATGCGTTCTATCCATTTGACACAATCACACAACGACAGCCCTGAACTTCCAATTCCATAAACCCTCGTGTGCGTTACCTTGTGATACGAGAAAAGTGAAGTCATCGTCAAGCCGCCCTCGATCCGCCAACATCGGTACCGAAGCACACCCACACCTCACGGGGTCGCGATTCTAACCCATCTCATCCAGCTAACTCTAAAAAGTTTAACAAATTATATAGACTTCCCGCTAGTAACAAGAAACCGTTGGAATTCATAGATTTAATATCCATGCCGATTGTAAATGACTACTCATAACCTCAAACAAACGTTTTTTAGCCCTCGCTAGACCGTTAAAACGCCTATGTCTACTAAATCCAAACATTAACTTAGACTACGCTAAAATTGCATTCCCACACAGTCACCTCCTGAGAAGGATTGCTAAGTTACGGAGAGTACGTGTCCCGTGTATGACTCCGCAACAGAGCCTAGTCACCCTGCTCCCAAATCACCCAAAATGAGACAAAGCATAATACAAATTGCTTTCAGATTGTCAATACCGCACCGAAGAATGCTACTCATAACCTATCTAAACTATTGAAACTTGGTGTTTATGAAATATGGCTCAGCGACCCCTAATGCTAAGGCCAAAATGTCTAAGACTTGTCAATCAAACGGCATGAAACCAATCCCCTCTTCAGAATGACCCAGTAGTTTAGGACAGATCCTTGCCCTTTCACATCATACCCTACCAACTAGAAACAGGAGGCAGAGTTGTATAGCAAACACATAAATTGTCACCCAACCATCCTTCTAACCCACTAGAGGATACTGTGCGGATTCGCAAATCCTCCTCCCCCGGTCTGTTAAACCCCTAACTTAGAACTCTTCCATTACTCAGCTCCTCAAATACGACCACTAAACCCCTTATCACGTCTGCAACTTACCGCAATGCCGACCATCCGGCATATCAAATTCAATAAAGACTTGCCCACACCAAACGCGTTAACCTAACACCCTCAACACCTTGCTAATTTCATTAACGACACAACAGAAGCCTCCATAGCCGCTAACGACCCCCGGGTAAAGTAGCTATCGACGTGGGAATAAACACTTGTACCTCAAAAAAGTCCAAATGCCTAACCATTTTCAGTACTTATACTGACAACAAATAGTACTCCCCCAACACCCCTTAAACTAGTTACTAACTTCACGCATTATAATCAGAAATTCTGGTCAACAGATCTCCATATCCTAATCCGTATGTAAAATGCATACCGAACCATGTTCTTCAGGGCACAGTAACATCAGAATAACCCTATATAACTAACTTTCCGTAAGTATGTGTGGTCATCTAAAACCATAACTTTCAATGCGTGTTACTAAAACACTTTAAATACGCAACGGATCCAGGTTACAAAAACACGGAAAGAACGATTCACTACACAATCACGACCCCTCCATGATCATTTTCCCACCCGCACATTCACGCACCCATATTGCAATAGTACCGTTATCAATGATTTCTCAAAACGACCTTCTTGATAGTCAAATTACAGGTCAACCTATTAAATAAACAAATACAATCCAAACACGCTTCACTCGCTTCAACCCAGCAAAATTACGAACCCCGGTCCAAGAATCTCATTATAAATTTCCACCTCTCCATGTTCAACTATTCTAAGTTTCAATGATGCCCTACACCTAACTTCCAAGACTAGAAACCCCAAACGAACCCCCTGATCCCATTTACCTCAGAATACAATTCGTACTCCGAAGCCGACAAAGACAGCCGTATAGCCCAGCTCGCAGAAATGACATCGGCGAGGGTAAATTGACAAGGGTTCTTACATTTCAACTTGACCCCACACCCGCCAATTCCTCGCAATTCGCGCACTCATATTCGTACACGTACCAACAGGGTTATATCCATTTGTTGAAAACCAGCTACTAAAAATAAACCCCTACACGCAGAGACTTCTATCGCCCACCACCTACATAACTGGAATAACAGGCAAATGCAAAACTCTCCGACTTTCGGCATTCGCTGACCCCCCAAATTACCCATAGAATTCTCACTCTGACGAACACGCCGCGCCACGTCCGCACATAACTTAGACGCAATAGGCAGACTCTAGCTCTTACATACAATAACCATTGTAATTAGATATGGCTCAACTGTTCTCGAAGGACAAATATCGCCGTTACCCTATCGCATACTCTACTCGATACCACTGGATTATGAAAAGTGCGCAACGCAATATGAATTAAACCGACGTATTCTCAAACTCACCCAGCCCACACTAATAGACAAGCGCAAATTTAACCTCTTATTCTCCCTTAAAGCTGGTTAATCCCAGCAATAGACACTCTAGTAAACTGAGAACGACCAAGACCCGTTCTGGAACAGAGATCCACCTTTCATACCGCCATATTCAACGCCCGATCACAATCGAGTAGTTTGCACTCGACGAGGTAGACTATCTTCCACCTCCTTGTTTAACACCTATTATCGATGCCGTACCTAAAGGTTCCATAAACTTCAGAAATCCAACAATGACTCTGTATTCCAAGCGCTGAAAACAAGATAGAGAATCATACCCAGTACATCACTCTAATACAGAGCATACGTCGAACCTATAGCCAGAACTTCGGAAACACACTCCTTCGAACCGTTTACCAACGGACGCGTTAACTAACTTCTTTCTAAAACCTACCAGCCTCTCCTATCCACTAATAACGATGAGTATAGATAACTAACGCGAAGATGGCAAACACACAGCCTGCGGCATTTAACTTGGCACACGACACCTGCAATCTATATATCAAGCTCCAACGTCAGTCACTGTCTCCGCTCTCAGTACTAATCCCACACAAACTATCACTTACGAAATAATCTTTGTCAAATACGTTGTACTGACTGGTCTACCTTAAAGTCTCTAACCTAAGTTTATCCCATCCTTGGCAAACTCGTTTTATACAATTAACCCCCACGAAGATGTGATTACGTGGCTTCGTATAATCTCACCGCAAATTTTCAAAAAAATCCTGTCCTAGGAAACTTATATCGTGAACATGCCCGAACTCTCTTCGAGAACCCCTATTACCAAAATTAAACACTAACAAAAGTTCATATGATTGTATTCAGACATCCCTTCTAAACCAAGTAAAAAACCCACTTCCATCAGGGAAACACATTTAAACTGTGACCTACAAACCTTCATACTTGTCGCTATACAAGAGACTAGCCATACACCAAACATCGGTCGACACCAACCAATCCCTACATGCCTTTTAACCCACCTAAACTGACGACAACCGCACACGTCGTCGACACTTTAACTTGTCAATCGCAATAAATACCGAGCTGCACACCGATATCAAAGAAACAGAAGAGCGACTCGGCGTCGCTGTCCCCCTAAGATTAAACCCTTTCTTCTTTAGCTTAGCTTCTACGACCCCATAATGAATTGTTAACCTACTTTTCATTACAGTTAACCGACAAATACTAAAGCGTAGACGACCCGCGCAAGTCACTCCCTACACAACTCTAATACATAACCTTTTAAGCTACCGTCTCACAAGTTGCAGTACCACCTGGATATCTAGCCGGTCTCCTCCGCACACAGCTCGCTCTAACTCATAGTGAATCAAAGTACCGTTGCCTCGCCGAACTCCACAACGGACGTTAAGTTTAGCTCCCCACACGCTGACCCCACCCGCATGGCGCCAAGCTTCAGCAAAAACGCATTCCACCCCCATGTTCAAAATAGTTAACCATTCCTAAGCATAAACATTCTAATGTGCCAAATGCATATTTCAGAAAACGACCCCCAGAATTAGATACAACAACATCACCTTTGAGCTTCCCATAAAACTAGTCTCGTCACTTAAGGAATAGGATCTTTAACAGATCTCCCACGTCCCCAACCACCCAATAATCTGCTCCAAAGAGCACCCGACCTACCTTAAAAATTCTGGTGACATTATACGTCTAAACGAATGGCCGGACACACGGTGACAGTATGCAGAATTAAATCCACACACTACCCTTCTCAAGCGTATAATAAAGTACCGTTTACACCCCCTTCCACAGAGATGGAGCTCTCCGAGCCATTTATAATTGCAGTAGCTTCATCTTTCCTAATGCTACACAGTATACGTCGACGATTACACCCAATAAGTCAATGGAACCTGCCCTCAACGCCTGCCTATATTCACTCCCTCTAGTCCAAGATGCACACCTAAAAGACCACCCCAATACAGTCACACTACTTAATACCACATACCCACCCATGATGTACGTTGCCCCAATTACCACACTCTTGACAAGTGGTGACAAGGAAGTATGTAGGCCATTAAATTTAAACTATGGTCACTCTACAGTACAATAGCATCCTGTTACCCGAGTGGCATATTACTAATCGTGTAGCCTACCCTGCCACTCTCGATCTGTATCTAACAAACATCCATATCACTTCCGTATACACGGTATCGTCCGCGCCAAACTCCAACTCATAAATTATGCCTAGAACACAATATTACTCTCTTTTTGGCAAATCGCCTCTTAATATCACCATCTCATCTTATAAGTCATTCATCAGATTATAGCTACGTACTCCTCGTATATCTTTCACTCCGCGACCCCAATGTAGTATACTAGCCTACAACAAACTACCTCTCTATCATTATCGCGAGGAATACTGTGACAGTTATAAAAACACATACTGAAACTTGGTCGTACCTTACCTTGGGCCTCAGACAAACAACAGGAACTGAAGTTGATCCACAGATTACCGGTCTCAATGTAGAAAGCCTACTACTGCTCAAATACCTTCAGTCGAAACAAA'

>> human_random = human(randperm(length(human)))

human_random =

    'GGCACATAGTGCAGTTACCGTCGACACCCCCATCCGAAAGACTGACCACAGATAGTGCACTCCCCCTAAACTCGCCGTTCAACAACATCTATATTATCCAAAAATTTAAAGTTAAAACATACCATAATTCTATTCACAGAATAAATCCTTTCACGAGACGTCCACACGCGTTTTCTCATCCCAGCCCTACATACCCGGTAAAATACCAACCTCTCAAATTCATAAGCGCACCCTCATCCCTCCAAGGGCTGCCTGCACGGCAAGCGTCCACTCCATTCCCTGGTCTCAAGAAGTCAGTCACTAATAACTCAATGATCGCAATCGTTAAACGCTGTCCGTAAACCCGCCTCCTAAAATCTTCCCATACGTCTCATAGCGCACCTCAACCTGCAATATCAGCACAAAGGAGCCGACGTCAATCAAAAAATCAACATGTTTCATAACTCTCCCTTCAAACTCATTCCTGTACACCTTCTATCACAAGCGCCCTCCACTTATCTACTCATAATCACTAACCATCGTAATTCCACTGGGCATGGCGCATATATTCGGATCGCCGATACACTCCACAAGTCCCATTAAACGCAAAGACTTAGGCTCAATTATAATAAGACCGTCTCCCCAAGGGAGATTCACAAATTTGCTTTACCTGATAGTTTCCCTCACCTAAATCACCCCCGAGCTGAAGCCATCTGAGTTTCAAAGTTAATTGCTAGTATATAAAGCTACTACAGTTCGAAAAAGTCCTTTANTCAACTTACAAAACTAAACCATAATATTAATGTACTGTTACCCAACATACATATCAGCGTTCCCCGCAGACCTATCTTAATCTACATCTCTCATAACCACAGGTCATACTTCAATGAGATCACAGTATAGTAAAACGTTTACTTCATCTCTAGTAAGATATATAATGCATCCACGAATGCCCATTGTACTACCAGACAACATAAATTCCGGAAGGGGATGATTGTCCCTTTCTACGCTTTCCTCCCGCTCTCAACCATCACAAATTAACTTAATTACCCCACAAGCGAAACTCGAGCCCGGAATACCGAAAGGCATCAAGACCAAGCGCGTAAGGACCTCGGCCGAATCTCGTTCAGTCAACATACCCTAATACAAACTCCCAAATACATCTGTCATGCCTCTAGTCACCGTCTGAATCCAATAATTTTCCAGACTGGCTCAACCCCCCCCTCCCCGTCAATAGTTCATAGGCCACTCCCCATCTGTCGTAATAAACCATCAAAGGAACAACATTCTCTTAGACCCTATATGAAATCACATAACCCCTCCGTACGCCGCCATAACACTATAACCCGACTTTGGTAGACTCCGCAATCAGTCTGTATATAATCAAATCCACACTCCCTCTGAACCAAACAGAACTATATTCACTAACATAATATAAATAAATATTAACAAGCCCTACAGTGTAGTCCACTACTTTCTCAGAAGGGTCTAACAGCTTACAGGACTAGACCCAATTATTCATAAACTCTCGCCCCGCAATTTCTTCTCGAGATAGTATAACGAGAAGAATAGAAACCACTATAAAAAGCTAATGAAACAACAAGTAACGTCCAAANACCTGTGACCCCACAACGCTGTTTTTTTCATCCCACCAAACTCTAAGTCGTCAGAAGTCTACCGTCAAAACAACCAGTCTAGCCCAACATTACTTATAAACCCGACAATTTCACACTGACACTCAGAAGCTTATACCGCCCCCACGTGTTCTAGCCTAAAGCTACCCCACCCACCCCTTAACCCATCAACGTACGACCAAACCCATATCCTAAACATTAGGAACGCTACCCGCAGAATGTCCCTATTTAGGAGCAACATTCAAGACAAATTGCCACCCGACTCACTCGTAGAAAACCAAAAGCACCATACCAAATCAGAAACTGTGTTTCCATGCATCTCAAACCCTATTTACAAGGGACTTCATAAGACCAAAGATCTTTCATCTCGACCCCAAAGCTTGCAAAATCTCGCCACCGCTCCCATATCCCAAATGGATTCCTTGTACAGACATACGACCGTTATCAACAATTTTTATGCTTTCTTAGCCTCATTCACACTTGATAACTGCATTCTCATACCTACACGCACTCCTACGGCGTTACCAAGCTCGAGCGGAACTCGATACCAAATCCTCATACCAGTGCTACTATGGCGATAAGCCCACACCGACCGTCTTTGAAACTAGGTACAACAATTTACTGTTCTCTGATAGCCCTTCCATAAGATTTACCTACACAACTGGTAAAAACGCACCACGGTCCATACCACGCGGTCCGTGCAAGGAAACCCTACATACAGACTCCGAATAAATTAAACTCTCGTACATTTACACGGTTCGACTGTAATGTTATACTACGCTATTCAACACTCATTAGAAACACTACCCAGCGAAGCTGAGAGACCTCTTCTCATATGCTACCCGGAATCTTTACGTACCATACATATTTTATCAGCTCCTTGCAAATGCAACACAACGAAGAATCCGTCCCAGTTCATGACTACTAAATAGTCCTCACTTTGCATAAGTTCAACTTATGATCTCCAAGTATTCTTCACCCGATTGAATATCAACCATCTCAAGTATAAACCGTTTAACCCTTAACCCCCACCACGGCCGACACAGATACTACAGAAGGCACATAAACACCTTATCGTTCATTTCACTGCCACAGAAACAAATATCGAATACAATTATGGACCACCGTAAAACAAAACAACATAAGCCAAATGTTCCACCATCATGAAAAAATACTAACAGCCTCACTACACAGGCCAGCAGGACTCCCCAATTCCTCATATAACTTAAAACTGCCATGCGTGACAAAATACCTCATCCACAAACCCACCAAATATAGATTCCGTATTTCCATTCAGTCAGAGTTACTCAAACCCAGCAAACACCAGCCCTAATTAAAAAAAATAGCCGACGCGCACATCCAACCCGGATGCTCCAGTGATCACGTAGTTTCTACACAAGTGGCGAACTTACGTAACCTACCCATACTAAAAAATAAGAAGCACAACCTCCTGTATCCGCAACAGTATTGTTTGCACTGGGCACTNAATACTTCGAAGATAACATACAATACTCGGACATTCTAGCTCCCCATCCGCGATTACTAATCGCTCTTCTATATTCTAACTGGAACCCGTCTGTACCCGCCTATTCACGGACTACCCCATACGTCGGGTATAGCATCAGCTATTCCTTCCAATAATTCCACCGCGCTCTCCTTAGGCCACCGCTGGTCAACACTCTTGCCCCACAACTTACAAGGACCATGGCCTAAAAGCCCGAAATGCGCGACTCCTAGCGTAGGAGGTACGAGTCGAAACTCAATCATCTAGCGCCCAGCTTTCTCATGAGCTCAAATGTTACGATTCACCTCACATTACAGATAGAACACACTGCCACCTTAGAAGCTTCCAACGCTCAGACAGCTATTCAGCTTTCAACTATGCTGCTGCACCATGTCTACCGGAAACCTTATACTCAAGCACGATTCCGTAATCCGCGTACAATCACCCCGACCATCCTCCCAGTTACTGCCCTGTCTACGGAACTGTTCTTTGCTAATAGCTCCCAATGGGTATTATGTCAATCTTCAGTAAACTTACTCACACTTGTTTACCTCCCCAAATGCAACAATCCTTCACAACACGCATCTCTCGACCAATCGTAGACCAAGTACCTACTCCCAGTGTTACCCTTGAACCTTCCATTGTGGTCCTGACCAAACACGCCACCCCACAACCCAACGCATGCCATAATTGAGACGCCCAATGACAATATCATACCATCCACTCAAATTACAGAAGTAGTCACAATCTAGCAAGCTGGGTTGGCAAGTGCCACCCTGCGCAAAAACACATGTAAGTTTGTATCCTATCCCACATTTCATCTGACCTCACACCAAATGTACTCAACATCTTTTAAGACATTCTGCCCCCCGTACAACCCTATCAAGCGGTATATACAGAATGCAAGCCTTCAATCTTATCTTCCGAACGCACAACAACCTTCCGTATCTACCACCAGTCACCAGTGTATAGAACCCGATTATAGGTACGCCTCACTACATTATGCTCCCCCATCAGCCGCCCCCAACGATTCTCACCCCACGTTTTAAACGTCTTAACAAAGCTTCCATCACCACCGTACTTATCTCCACCTAATACGGCCCCTCAGACACGACGCCAGACACTCTCACGGGCAATTGTCTCAGCATTATGCTCTCTTATACCGTCTCTCACATCGAATTAATCTAGTTCACCACTCTCTTCGGATCTTGAAAATCGACTAAGTAAGAACACAGACTTAACCGAAATCGCCTTGGTGCAGCAGCTTTCCTNCGCAGAAATTGCTCTTCCGTAGCACTAACGAAGTCCTTTTCTCCATCTCCATGCATCTCTGACAAACAAAACACCATGACCCCCCAAATACCCCCCTAAAACTGATTCCCACCCCACACAACATACTGAATCGGGAAACTCACACCACGCCTAATATCCGCCTAAGCATTCTGTTGTACTCAGTAGGACTGTGTTAGCAACAATAGAAGCGACAATATACGACCCCGCTCGAATTGCCATACCTTACCTAACCACCTAGATCACCAATCCAGCCTACCACCAAGCTCGTACTCCCGAATACTAGTTTCCTATAGGGTGAAATATGTCCCACCATCAATACAGCTCGAATTTGTGCCTGTCGTCACACTTGCTATGCTTTATGCTTAACCATCAGTTCGCAGACTTTCACGACATAACTTCACTTACAGAGAAAAACGTCCCTGAACAGTCTCGATAAGCACAAGGACACCCGCCGCCCTCAAAACCCCCGTCGCATTATAGTCCCCTGCTTCTTATACACATGGCCCCCCTCCACCTTCATAGTCTAACNCCGCCCAAACATCTTCTGTCTTACCCTGGTACGACAATCCTAGTTACTGAAAACTCCACCCCTACATAAATTGCATCCCCGCGCCCACAAACACGCTCCTATACAAAGCAAAATACACGATGGGCCTATTAGCTCTCTACTCTCTTCAACAAAGACCGATTACCCCAATCAAGCCTCCACTCCTGGAGGCCCAAAGACTCTTGAAATCTGATAGTATAAATTACCAATTCTACAACCTTTTTCAATTTGACCCTCATCCATCTAAAGCTATCCACGTTCAGTGACACCCGCCCAGCTTCCTGATACCCTTACCGCATCTTGTTCACGTAGGTTCCTTATAGGTTCCGGACTGCAACGCTTTTTTCTCATCCATGCACGCCACCCATCGGCCACTGACAATTGTGATACGCTCCATATACATCAATACTAAACGACACCTCCGTAAACCTACCCCTAACAAAATTCCATCCACACAATAATCGTCTCGCTATTTCAACCCCTTCACCCCGCACGTTAGACCCTCCCCTGACCTCCTATATTAAGGAAGCACTACTACAATGCCTAGACTCCCCCCTGTCACCTTCTACTCGCCAATCACCACTCGTAGCATACAGAATCCTCGTTTAACTATTGCTACACCAAGTCCGACATCTACACCGATGGCTTTCTCCCCCATTTTAGCAATCCGTACCGCCGATAACAAACTACGTCACCAGTAGATCAGACACGGGAAACCATCCCATCTATAAGTAACCATCCCACCCCGATTAAATAATGCGTTCATACAATGCATCCCAGCGACGACCTGAACGTCTACGTTACTACGCTTAAATGCACTTCAGTGCCCGGCAACCACTAACACGTGCATCGCACTATTGTTTTCTGACGCAATTACAGAACGCACTTAGACGTTTTATATCACTACTACCTATTGCTTACTGCATTAACCTCTTCTTCCGTTCGGTAAATCACACGTTCCATCAAGCAAAGAAATCTAACTAACACTCCATGCCGGTACCCTTCTTTGCCATGTCCTTACTCCCCAATGCACCATCTCACCCTTCCCTCCTTAAAAACCGTTGAATTTAAGCTTCCCCCGAACTCCCTAGTCATCCCCGCCAGCAGCACTGAAAATTACTATTAAACACATAACCGCGAAACGATACAATCCGCCAATGTAGCAAACTAGGAACACCTACCCGTTGATGTCTAACCAACCAGCTAGTTGCGATTTCTTCCCCCGAAGAACCTAGAAAGCGAACATTTATTTTACAGAAAGTACTTTCAGCAGCCCAAGAGTACCCTTCTCGTACGAAATCTACAAATAAAATTTCTTGTACCTAAGATTTACCCTTCTGTAGAATAACACCACCCCATAAACCTTTTCCAGAATCGGTATATTTTTAATATTTGCCTATGCCCAAGTCTCCACAAGAATCTTCACGAAGCCTCTACCTATAAGTAGTTATATTTACCGACCTATTCCCCCAATACCATCCTATACTACTTCAATTCAAAGCAGTACTTGACCAGTACCCGACAGTCTTGCCGTGAATCATCGATGCTCTATACCCTCACAGTCAGAACCATACTGGCCATCAGACCCGATAAATCGCATCATCCGCGTAACCCACCAATGGCAATGTACATATCTCCGTTCGATCCTCACCTAACGAAATTATTTCCATTGAATAGAGAGGATATATGTCCATACAACCCCTATGAAAGATTACAGACGCATGTTGTACCAACTTACCTGCTTGCACAACCACACTCAAACCCAGTCTCATGTTCCACAAGTTTACCTCAGACCTTGCTATAGGATTCATTGCCAAACAACGCACATACTATATTCTTTTTCAAACTCTTATCCATTCACTTCATATTACTCCAATGNCTCCTTTGACAGACAGGTCCACTTACATTCTTATCAACCCTGTTGTAGCTGACGAACTAGGCTTACTCCGAAATCTACCCAATCATATAACCCTTCCCCTGGTGTTACATACTTAACCTACAGAAACAAATCTTCAACGCGTTGGGCTACCCGCCAATTTCTGCGAACGTACTGCATAAAACCTTACCGCGGAACAATTACCGACACAATTTAACAGCCGGCCTCCGTTAATACAACCTTATCAATAGTAAGCGTTACCAAGATTCTCTTTCCTCTAACCCAAACACACTTTGCTTAGTAATCTTTTCGCCCACCAATTATACCTTATTCCCAAACAATGGGAACAATAAGCACCTAGCTCAGTTAACTAAACCAGCTTAACTAACCTCTTCAGTACTGCTAAACCCTACGACAAACCATACTTGCACAACCCACTGACCGTTCTCATCGATCAGTTTCTGTCCTCCGAATACCANCATTTCTGGACATTACCTACATGTCATAGGGAACGGGAACAAGCGCACTTACCCCCATCTAATCCCGCCCATTACAACAGAAGAATTAACTCCCACGCTTGAAGTCATCCGTGCACTCATAGTCCCACGCAACCCATGGAAACCCCGAACCGCTACAACTTGCCCTGTTCGTACAGTTCCATTACTACTATCCTCTGACTATAATACAATCGTTTCTCACCTTTGATCCTTCTGTATCTTGCAGACACCCCCCTCAAAACGAAACAAATCTGGAGCATATACCTGAAACACCGGACATCTTAAAGCCTCGCAAAACCAAAACACTCAATTTATTCATCCACATTCCTGATAGAAACATTAGAATGAGACTCAACATCGCCAGCTTCGTGGCTACGAACCCACCCCGCACTCAAAAACTAACTTATTGAATTCAATCTACATCTACTTTTGAGGCTACTCGATCGGCCAGATAAATTATTACGCACCATTCGCACTAAACCAGCCTACACCCCAATTCGCTATTTAAAAAGTTCAATTGCACCTCCGTGATTTTCACCTTCCCCCTATCAGTCAATTCCACATGAGTCTCAAATCCCCGATCTACACGCAACCTATTACCCCCTACTACTGCATTTGTTAAGCCGTTGTATCCAACAGTAAATTTACCGTAACTACTTCCTCTGCCCAAGCAGCCCTACAAAAACTTATACACCACCCCTTTAGTCGTGTAACTCGGTCTCCACCCGTCAAAACTACATACAAATGACGAGTAGGTAGCTATATCATAGTATAACATTCCTTACATATACGTTTCAACAAGTCCACATCCTTATCACAAAGTACTATAAATTCCAGAAACTAAAGGAACTCACGATCCCGTAGTCTGAGTAACAAACCGACTTTTCACGACCGTACTAGCCATTCTATCTGTACATGCTCCATATGTGTCTTTGCCTGGAACAGCCAATTGCTCCCCCCGATTCGTCATAGCGCAACACAACAATACTACGATGTACCCCAACTGCTCACGGTCCTATCAACACTATTTAGTGCATCTACCGCTTCTACAAAGACTAACATCAACAACATACCCAATACTTACTCTCAAAATAGTCCGCTTCCATACTACTCAACCAAAACATACTGACTTTAACCTCCCCAACGCTAAAAAGCCCCCAATAGACCACATTTTCTCTCCCAAATCGTCCNAACCGAATTGAACACATTTCATCAAACCTCTACCGGGATCACAAATACGGTATTAACCGCCGGTCCCTCAACCCTCTACCCTAATAGAACAATAATCGCCCTAAGGACCTCATAGTTCAACCACATTGGTCCACACTACTCACATTAATTTCGCACCTTCTACAAGTCTTGAATCTGAGCTCGACCTCAGACTTTATGCAGTGTCCTCGACTGCGTTTATCACGTCACCTCTGCACCCAGGGGAACCGCTACACCTATTCGGACCCCCAAACCTCAAGCTCCAACTCTCACACCGTCACTGTGCGATAGCGCTAAGCACACGACAATACAGTACTGTCTCCTATATGTGATCATTATACTAGCCAACAAAAAATCCACCCTACATTTCGACTGATCATATCATTATACCCTATACTCCGTCAACTAAAGTAACACTGCTGTTCGACAACACGAATTATACAGTCCGTTTCTATCACCAATATCCCAGACCACTGCGCCAAAACACAGGATACTGACCTTGGCAAAAGGAGTTTCCCAGCTACGAACTAACACAAGTAAGCAGCGATGGCCTCTTCACACTCCAACAACGCAACTACTCCATAGCCCCCATTCCGTACATCTAAAGATTACCAGGACCCCACTAAACCATAGTCCCATTCACAGGCAAGCCACCTGCTTCCGTTCCCCCAAACTCGCCCCTTAGTACCTGCGCCTATTATTCCTGTTGAGTCGTCTCCTTATAAGCTGCTGTCTAGGTTAAAAGGCACAGAAAACAAATTCACGAGACACACTGCTACAATTACCTGATAATACAAACTTACCTTTCATGCTCCCCCTCGGTCTTCAGAACTAAATCTAGTACACCCTGTCGCAAATTCGCTTAGCCAGCCGCTTAGCACGGCTTACGGCTCCCACATTTATATATATTTTAAAAATCCCTCCAAAAATATAAAAGGACGTTTATTGAGAATAGCCCCGTCACCCTTAGTGGAAACAGTACTTAACCAACTCACGTGCTAGTACCTACTTATGTATTGTCACCTTAACTAATAACCTAGCATACAACGTCACTACGACCCTCCTAGATTCACCCACGAGCACTCGCAACGCTCATAGCCCAAAGCCAAACAACAAAGCCTCTCCGGTGTAAAGTGCATCGAATCACATAAACAATGACCCAGCCTATACACTACACCGCCCGTCCATAACCTCTCGCCCACTGCCTTACTACTAAGACAACCACACCAGCACCGGGGTACGAGAGCGATACCTGACGTCCAAACAATATAAACCATAATAGGTTCGCTGGCCCTATGCCATACTGCACTACTACAATGTCTCTCAAGAAATTAATCAGTTCCAAGCAGATTTGCNCCAGTTACTACCACTACCTGCCACGAATAGTAGAGAGGATACAAACATTTCCGCAAGACTGAAATGGCCAAAGCCTTCAGTCCCCATACGATATCTCTGCTGTACACACGAAATGTGGAACAATCATTTACCTCACAGCGAAATTCTCATACCAGCAAAACAACCAGCCTTTCCAAGAACCTACAGCTACCGGCTCACTAAGAACATTGTGGAGCTTCAAGAGACAGCAGAGTCAATTGTAAGTGAATAAGCTCTACGCAATGTGTAAAAATTCAACGTCTATCCCAATTGTCGTACATACATATTTTTCTATACATTTCAAGAAATCACCCCACGTCCGTCGCGCCCCGGATCGTAACCGCCCATAAATCAAGTTCAAGGTCTACACAAGCTCCACGAGATAATCATGGAGGTAAATTTCGCACGACTCAATTAAAAACATACCTTTTGCCTGTCGTTTTTCTTCAAGTATAACCAATTTAATGTAAATCACCCAATTTGGCACCACACGCCTCCCCCACGATAGACAATGCCGCNTCATGATATTATCCTGCACTCATCACATATTTTTTGATACTCAAAATCCTAAACATCGTCCCACCCTTTATTATTGTCACAAATGTCCACCCCCGTATGCTCCTTCCACTTCATAAACCTACACGGCGGTCTGCTTTAACGCAATCACACAATGAGCCGATAATAACCATTCTCATTCACCCATTTCTGCTGGATCAACCCATCTTTCGCAATTACCCCCCGATTTGTGTGTGCATTGATTCTCCCAGCCATCCTATCTCAACCATCCCCCCTCTCGTATGCAAAAACCTTCTTCACGATTGCACACTACAGATTTATCCTAACCAAGGCGCCGCGATCACACAACAATCGCTCACCCGAGTTTCCAAACTCTACCCCATCCGACCGCATAAGCTTACTCACTTTACGCAATCTCAAGCACTTCCATGTACTGCACTAAACTTCCCTCACACGAATAATCCACTACCGTAACAGGTCATTACAAGAAAGTCATAACCAGAGCTGTTAAATTATACTCTGACCCTTACCTAAAACTATTACCTATAAGTCACCCTACCCGGATCAAGTTGGTTTACATAAATACCAACGCATAAACACTGAACAAGATACAAAAGTATCTAAGTATAGTGAATTACAAGCATAGATCACAAATTACGACACCCTACCCTCCAAATCCGAACCCCCGCGGACTGAAGAGATACAGTTTTTCATTACCCTGTGCACCCTACGCCTGAACAGCCTGATATCAACACATTATGTCCTTCTCCTATGTCTTTCCCACACCAACTAGACCCCAACCACCGGAAAGGGAATCGGTGAATAATCACTAACACTACACAACAATGCTGTTACATTAATATCACTTCTATCCACATTAACAATCGCCACTCAACGCCTTATCAAACGTACAAGACAACAAACAAGAACAAGTACAATACATTATCCCGCTCGACGAAGACGCCCCGTCAGAGTCCATAAAGTTAGACCACCAAACTTAAAACAACGGCAACCCTGCCAAAGCCCACAGATTACGCAATTACACTACACACGTGTCACTTACTACAGACTGTTAAAAGGCAAACCCCACTTAATCACTCAAAAAGGTGTATATACTCCACATACTAAGGCTTAACATCCATCTTTATCTCTACCAAAAGCATAAACGAGCTCCAAAACACTAGACCTGAAAACGGAACACATCGGAGTACTAACCACGTCTAATATTTGTAACTCACCCGCATTCAAAACATACACATCCAACCCGATATCTTAAAAAAGCCCCTATCCGCTCATCTCATCTTTCGTTTCTTTAAATGTCTCCCCCTGGATCAATATACTACAAATNATAAGATATCCATCAAACAGGGAACACTAGATCCTCCCAGCTATTCACCCACTTCCCTGATATTAGCCGATATCCTCATCTCGCTTGATCTTAGGCGTAGTAACGTTCCCAACAGATAATTCAACACACCATGATTTTCCCAAACTTAAAATCAATACTTCCATGGAACACCCTAGAACCTAGTAGATACTGGACTTCCATTAGTAAACTAGCACCACTTTCCACCCACTAACTCTAAAGCCGTGATTTTTCTAGAAATACCCCATGGAATTCCGGTTGTACCACCCCAAGACCACAGCAATACGTCATCAATCCAAAACGCCACCCCTAATACTTCATGCATCGCCACTACATAATGCATTACGATTGATAACATACATCAGGAACCTCTGGACCTAAGACTTCGCCTTACTTAGCAGAGCGCAAATGTCCCCGTGTTTTTCTCTAAAACTCCCTTAACATCCGCACGTGTTATAGTACTTTCCGACTAATACATCGATTATTCTACCAACGTCGACAAAGGCGGCCTTTCCGCACGGTTCTTTCATGTTATCTATTGTTCAACGCCCCTTCAATCCTTAATAACATGTCGTAACTACTGATTTCATACAACAAGTGGCCAAAGAGTAATACCCTCTGCCTAGTGAGATATGAGCCAATACACTTGAACGGGTATATCAGTTTCCAAAGCCACAGCCGTTTACTCGGACTCTTCTCCTCGTGCTCCGCCTTGCTGACCTGCGCTACCAACCCCCTCCTCTTCACACTATTATCATTATACTGCAAATAACAAGCTAGTGACCGACGTTTTACCGCAACTAACGTCCCAATTCTCTCAACTTGATAATCTCTCTGATATAAATTATGTTTACGTCCGACCGACAACACCCAANCTGGTAATCGAACTTTCTACTAACCTGACACTAAATTGACTTACATGTCCAGTATCAGCCCTTATACCTTCTCATGCTAATAATCCTTCCAGTCAGATCGCCCCTACCTCATTGATGACCTTACCTACGATCTGGTAGCAAATTTTAAATCCAGCCCCCGGTCACATTATTGCAACACCCCGACCGATGAGAGAATGTAACGTTAACACATCTAGAAACGTTCNCCCCACCACTATCGGTTACAGGAAAGATAATTATGAGCCCACATACAATAANCCCGGAAACAAGAATGCACTCTCCTCAAAAGACCGCCTGAACTCCACATTAGTTGCATCAGACCAGTCCAGCCTTTGCATACGCGATCTGCCTGCCACGAACTAATTCAGATAACATAAATCTATTAAAAAATTACCTACTATTAGAACCACTCCGCCCTTTGCTCACTCGCAGATCCCCCCTTTTTTAAATGACATACACCTCCATCTTTACGTCCGTACAGTAAAATTCCCCTCCCAAAATGCGACCAGGCCACACCAGAACATAATTCACACATTTATGCGCTCTGACTACACTTCTTCACAACCAAGGACCTGTCTTAGTCCAAATTCATGTCCTCACAACTTCTTATATAGAAATCAAACTCCTCCAACTTAACATTCGCCAAACAATTACAAACAGACCCGTGCCTTACGTCTGATAAATATTACTCACAGCNCAAACGCTATCAAACACATTCACGCAAAGCCGACCCTATGAAAGTGTTAACCGAACGCCCGACCCAGATCAGACATTCCAACGTAAAACTTAAACCTTCCTATCGCCTTCCCCATAGCCCGTATAGTTACTAGAGATTAGTCTCGTCACATCATGACTTAATGGCATCCTGCACCCCCACATTCCACCATGAATTACCTACACGTCTCATTAGTGTAATAACATGCCTCGTCATAACATTTTTGATAGCAAATTTTCAATTTGTATGGCCCACGCATCCATCAAAGATAATTATTATTTCTCATGACTGTATCTTCATAAATCTGTTCACCCTGCCACGGCCTTATCTCGTATATCCGAAACCACTGCAGCTTGGCCATACTCACTTTATGCAACCCATCCTGGGACACGCCCTACGCGTGGCATGCTCCATCAAACGTGCGGCATACCTACCTCACCTCCACACCTCCCTCGCCCGGGACATACTTGGGTAAAACCCCCGCGTTTGCCCCCCACAACACACCTGCCCAAAAACTCAATTACTACTGAAACTCAGGATGTCGCCCACTTGATAACACAGTATCTTAACACAGTCACCAGAGATAAATCCATAAATATTCCCTCCTCAACCTTCGAATACCGGTAAGCTTCTATTCNCCCGCAGTAGCAAGTTATATACCAACGGCACATTGCGATCCCCCAAAAAACGCCTACTCATATGTGTGTCCGCGTCCCCAATAAAACCACATGCTAAAAACCTTCTGTCCGAATATAAGTATCTNGGTGTCTTTTATCGCGAGCTAACACGCCCTGAAGAGTAACAGACACCAGAAAACTATCAACACAGGCCCAACCCGTAGCCTACCTGCTCAACCCGAAAATACATTTCAAATCTAGGTCTGCTCCGTAAGCTTCCGACGGGCGCGTTTGCACGTTGCTAACCCTACGCCTGAAATTCGCTATCGCACACCGCGTAAAACAACACTACGGGCCTCGCCGCAGCGCTTACGGAGATCTAACGGTCAGGCTACGCCTCTATTACACCCCCACTTCTCACATACGATCACTTCCGCCTTAGGCACTTACACCCTCTAATTTATCACATGAAATACAGCCATAAGGCCACACTCTTTTCCCTAGCGATCTTCAGCTTCCAGCGCCCACCACTCTAACTGCGTCTCGATGTTCCTATTACAGTTATATCTTCTACCTTCAATATAACTGTTCTGAGAACAGAACTCAATGCTCAGAACCACGAGTCGATCAGTATACTATAATAAATCTAAATTAAATTCCCATGCTGACCCACTTAGGACGTCCCAATCCATGCAATCAACTTCCTCAACCCGACAACTTTCCATTCTGAGCGCCGCCCGTCCGCTACATAAATCCCATAGTTTACACATACTCCTAACTACTGATCACATCTTTATGTCACCAATCAGAAAGACACGTTCGCATTACTATAAAGATTCTGCACTGTTCCCTCCGTGTTCCAACCTGTCACGCATGGAAACACTCTCATACCAAGGCTTAACTAGTATTCATTGTAAAGATAATGGCCCAATGCTTTCTCCTTCAAATTTCGTCTACGTTGCTGATTCCACCCTCAAACCCGACTTACCATGCTTCTATTACACAACGAACCTTAAATCGGACCAAACGCCTCATTCAGTCCGCGACATCTAAGTTCCATATACTTAACTTGTTGAGACCACAAAATAGAACCAGGTCCCCTTTAAATTAACAACCCCTATCCACCGCTCGATGCCTAATAGTGCCTTAAATGTTCTTAGCACCCTATCAGAACTAAACGATCTAATTAAAGTATCATCCCAGAACCAACGAACATTACCCCTCACCCATACACCAGAACATTCCAGATCCAGAAAACACTCATGAGCGTTAGATGCCCCACCATAATCGCTAACCCACACGAAATAAATCCTGTCACATATCATATCAAAGCCACCTAAGCCCTAGCGTCCTACCCGTACCCCGGACAAAAGCCAAAATCACACTGGCCCTGCTTTCCGAAATTGTCAGTATGAAAATTTCACCAGTATCACGTTTTGCAGATGCCCGAATATCCCCCAAATTAGTCCAACCGCCTCCCTGAAACCAAGCAACTACAGTGCATATATCGGTTTAATCACCCAGACGCTACTACGTAATAGCCAAGACCTCTCCTTGCTTCGACCCACTAACCCACTCTTCCTTTTAGCGAATCGACAATTTCCTAGCACCAACCCCCTCACAATCTCACCTGAACAACTACCAACCTTGGACGTAATTAAAAACTTTCGTATCCATACAGGCCAATTGAAAAACAATGCCTCGCACCCGTCACTAGCCAACACTCCATCCCGCAGGGCCGTAGGCATCACTCACATACTTTTTCCGTGCATTTGCTTTACAGGGCAATCCGCATTCAATTAGCGAGCCCATCAATATTAATAGAATATTACGGCTCTGCTATTCCTTCATAGCTCCCCCAGCTAGCTACTCTCCCCGACCATTAAACAGACCAGGCTGAAACTTTTGAAACTAAACCACAGTACGAGCTACAGAAGTGACTCTGAATCGGCCAGTGGACTACTTGTTAGTGCGCCTTCCGGACAACCTACTACCCACACATTCCATAGGAGAGGCTTATACTCACCTATATATGGAATCCAACTCCCACCTCCCTAAACACCAACGAAATATACTCCCACACCATGCACGGAATGATGTTAGTCTGGAATCCCTTCCTCTACTATCTACAAAACAGTGCAAACTAACAAGTCCCCCCCCAACGCCTAAATTCAAAAACAACAACAGCCGTTACAAGGGTAATTAACCCAGTTAAAATCCGCACTCGAGACAGCCCATTATCATGCTAGACTTCCGCACCAGACAGGACAATCTCGCATACTCCTGTGGGAACTTATCTCTATAGGGGTCCCATTAGCTCTCTCACCCCCAACATTAATATCTGAGTTAGCATAAACACTGGTCCTACCCATAAAACCGGACCCCGCNGAGACCAACGTTCCCAAACAATAAACTTCTATAAAGTCCATTTGCATGGAGTTCTTCCAACGCGTTTCCATTTGTAACCAAACTACAAGTTAAAACAGCAATGCAAATATGCATGACTCCCTTCAATGCACTTGAACGTTGTACCACTACCAATTCCCAACCAC'

>> [score, localAlignment] = swalign(hippoPORF, humanPORF)

score =

   22.3333


localAlignment =

  3×157 char array

    'ILISLILFIGSTNLLGLLPHS-FTPTTQLSINLGIAIPL-*AGTVIIGFRNKTKISLAHFLPQGTPTPLIPM-LVIIETISLFIQPI-ALAV-RLTANITAGHLLMHLIGG-ATLA-LINISITTALITFIILVLLTALEFAVAI-IQAYVFTLLVS'
    '||  : |  |  |::|  |:: : | :: :|:|    ||  | ::|  : :   ::|:  |   || | ||  || ::   |||    :|||  :  :  |::  : |||:  ::|  |:  :| |:| :  |::  :::::: |  | ::  || |'
    'ILGYIQLRKG-PNVVG--PYGLLQPFAD-AIKLFTKEPLKPATSTITLYITAPTLALTIALLL*TPLP-IPNPLVNLNLGLLFILATSSLAVYSIL*SG*ASNSNYALIGALRAVAQTISYEVTLAIILLSTLLISGSFNLSTLITTQEHL*LLLPS'

>> showalignment(localAlignment)
>> [score, globalAlignment] = swalign(hippoPORF, humanProtein)

score =

  359.3333


globalAlignment =

  3×226 char array

    'MNENLFASFITPTILGLPLVTLIIMFPSILFPAPTRLITNRLVSIQQ*LIQLVSKQIMNIHNHKGQT*TLILISLILFIGSTNLLGLLPHSFTPTTQLSINLGIAIPL*AGTVIIGFRNKTKISLAHFLPQGTPTPLIPMLVIIETISLFIQPIALAVRLTANITAGHLLMHLIGGATLALINISITTALITFIILVLLTALEFAVAIIQAYVFTLLVSLYLHDNT'
    '||||||||||:||||||| ::|||:|| :|:|:   ||:|||:: |||||:|:|||:::||| ||:||:|||:|||:||::|||||||||||||||||||||::|||||||:||||||:| | :|||||||||||||||:|||||||||:|||||||||||||||||||||||||::|||: :|:: ::|| | ||:||| ||:|||:||||||||||||||||||'
    'MNENLFASFIAPTILGLPAAVLIILFPPLLIPTSKYLINNRLITTQQ*LIKLTSKQMITIHNTKGRT*SLILVSLIIFIATTNLLGLLPHSFTPTTQLSINLAMAIPL*AGAVIIGFRSKIKNALAHFLPQGTPTPLIPILVIIETISLLIQPIALAVRLTANITAGHLLMHLIGSTTLAISTINLPSTLIIFTILILLTILEIAVALIQAYVFTLLVSLYLHDNT'

>> showalignment(globalAlignment)
>> [score, localAlignment] = swalign(hippoPORF, humanPORF)

score =

   22.3333


localAlignment =

  3×157 char array

    'ILISLILFIGSTNLLGLLPHS-FTPTTQLSINLGIAIPL-*AGTVIIGFRNKTKISLAHFLPQGTPTPLIPM-LVIIETISLFIQPI-ALAV-RLTANITAGHLLMHLIGG-ATLA-LINISITTALITFIILVLLTALEFAVAI-IQAYVFTLLVS'
    '||  : |  |  |::|  |:: : | :: :|:|    ||  | ::|  : :   ::|:  |   || | ||  || ::   |||    :|||  :  :  |::  : |||:  ::|  |:  :| |:| :  |::  :::::: |  | ::  || |'
    'ILGYIQLRKG-PNVVG--PYGLLQPFAD-AIKLFTKEPLKPATSTITLYITAPTLALTIALLL*TPLP-IPNPLVNLNLGLLFILATSSLAVYSIL*SG*ASNSNYALIGALRAVAQTISYEVTLAIILLSTLLISGSFNLSTLITTQEHL*LLLPS'

>> showalignment(localAlignment)
>> [rscore, rglobalAlignment] = nwalign(hippo_random, human_random0
 [rscore, rglobalAlignment] = nwalign(hippo_random, human_random0
                                                                 ?
Error: Expression or statement is incorrect--possibly unbalanced (, {, or [.
 
Did you mean:
>> [rscore, rglobalAlignment] = nwalign(hippo_random, human_random)

rscore =

   1.7085e+04


rglobalAlignment =

  3×17571 char array

    'CGGTAACCGAAATGCA--AACC-TAGACCCAACCTCACAAATTCCGAGAATC-GA-TACTGTTACCTATTCATCCAGTTTTCTCCGCACTTTCGACTGTAT-AATCACACCAGCGGTACCCTACCGGACGCTCCTCGACAAATTATCTACTACA-AGACTATTAATCC-CAC-CGATTCATCCTTCTAAAAAACTGTGCTC-TCTC-GCCC-A-AT-CTAAAAACCTTACCTAACTTTTGAATA-ATATTCG-ACAGGTCGATAAAAAAAATGCTGATTGACTATTCTATACAATAGGACATCTCCC-AGCCGCGAGTA--CTATCACTGATATCACCAAGTATCTCCTACTACTAAAACAATAC-TATACATCGGCTCACTAGAACAGAATAATTTAAAAAAATACTCAAGATGTACCGGCTCATCCGACTTGTATCC-CAACCTTTTTACCATTTACTAATGTCAATGGTGGACGAACCTCATTGACGAATGCATCTAGCGGACTCATCTCCAATAGCGCC-ACAAACTCCATATGGTCCCTTCTCACACATTACACCACAGTAAACGATTCCCATTAGTAAACAGCCGATAGTGCTAATTCACCTTCTACTTCCCCTACATAATTATTATCTCCCAACCGGTACCTAATTGCTTCTGGCACATGAGTTCTCGATTACATCACATAATTGTCCCTATATGCGTGTCAAGATCATGCCGGAGATTAGTATGGGGATCGATTAACATCAATCAAACCTCAGGAATCAAAACTCCCTACACTAACTTAGAGACCTTAGCTTTTTAGTGACACAA-TGCAAAACAACTAATATACCGAATATCCGACCCCACTACCGCTACAATCTAGACATTAATGCACATACTAGTATCCTACAAATG-AACAGCAAACCCATTGCTACAACGCTCACAAGAACTGTTTCTATAATGTCCGCCTTCCATCTTTCGTGAGGGAACACCATATCTGTAGAAAG-ATA-TA-AACGCAAAATAG-T-CGTTTTAC--CCAGAGATACTCA-CG-ACTACGGAGGCCAACGCAACCCACGCCTTCAG-CGATTGATGCCACTATACGGGTTGCATCGGTCACTTAA-CTTCAC-ACTT--CTAAAGTTCCGCCGCAA-CG--CCTC-TTCTCC--ATTGACTAAAGGTCCCCCA-CCGCAA-CTCAATAACGTACACCTATCCTATATCCGTAGAACAAACAGATCGT-CCCTCTTTCATTACCAACCCCAATTGAACTAACA-ACATATAG-C-CC-T-AG-GCCCACT-A-TAT-C-G-CT-GCA-ATACGTCCCTAAGCAGGC--TAGTTCCAATAATTCAATCTAACTTAATGGCGGTA-ACACAAAGCACA-ATTCGAAAGTCAACAGAG-CCCCGCA--AACACTAATAATAGAACAGTATTTAATCACATATC-AAGACCC-CCCTAGCTC-ACTGCCGCACTAAAATAACGGATGATAATAATTGCCCCAGAGGCGCTTACTCCTAAACTAAGCCGTAAGAAAAAGTACTTAACTCCTGTAAGTCAAATCACGGCC-GCAAATACATTCC-CT-C-CCCCCTAACCGGTCTTGTCGCCTTTATAACGTAACTGCAA-AAATCAT-CGCAC-CGGCTCCCTGTCTG-GTC-TCTTGAGCACGTTACCTCACCAGGACCTGATTTCAAT-TCACTGACTGCC-TATCCGCCTCT-TC-TCCAAGTA-GTATCTCCTAAAAAAATAGAGACGGTC--CCCACCGCACT-ACCGTCCACTGACCTC-ACTG-CCTAA-ATCC--T-AAGCACCAA-AAAACTCGTTTACCAAGCTA-TGCTAGCT-ACCCTTTGTTGTACGAAAACCTATCTGTCGATG-TCTAAACGACTGATTAAAAACAACATAAACACTCATCAC-T-TTCGTAGATACAAAACGGATCTACCTAATCAAAGATATACGGCATATGCG-GGTATG-CGCTTGTCATGTAACAAC-TGCCAG-C-GCTT--CAAGCGTCT-GCGCACAGAGGG-CTTCT-CATAATTACGAAGCAAGTTACC-TATGACCCTAC-TAGGTAATCGAAGGGACCACATTCTTCCCCA-CCGCTAA-AGAAAACATGTAATAACCGAATCAT-CAAACACAAACAGACCTACCAAT-CCGTAACCCTCA-AGACTACAT-AAC---CAACCG-CTTAAAACAAATTAGTATATAAGCTAGAAAGGTGACCTGCATAAACTCGTCATCTCAGCAATCCCGAATAAAAAC-CG-T-CTAC--CATTGTACAAGCTAGCAGTCAGAT-TCGAGACAAGCATCCTCAAA-CTGTATCATAATTACAGATTAGCGCACA-AGAGCACAAACGAAACTC-GT-C-TC-CTAAGCGGTGC-CTTAAAGAATCAGC-TAA-A-TGACC-A-GCTAC---CAAACCCACACCTC-CTCCATGCAAATAAGT-CATAAACTATAGCCTGTATCGTGATGAACGTACTTACCGGTACAAACGATTTCGTACATACGTAACAGTA-AATTCCATTCCCTTAATAC-CTAAACCATTCTCAATACAACCTCTT-GCGGCAAACCTAAAGTGTCCACCTC-CAT-T-C-GCCTTG--TC-TCACCAACCTTA-ATGGTTCCA-ATCAACAATCCTAACCTCCACTAAAGTA-AATCTCAATCTTCCAGATCTACGCGATATGATAGGCAACGAATCTTCACACATTGGACATTCGCGCACTCTTATATT-GCCTACCGA-CGTATTT-AA-AGTC-CATGGAT-ATCCACTTAAAACACCAATCCCTGCTC--CCTAATCAGGCGCTAGAAAATAAC-CTCATAGACC-TA--G-ACA-CTC-CTATCACAAATAC--GCAT-G-A-AC-GTTAT-TACCTGACCGCATTCCTGACCGG-A-AATCCATGATTTATACGTTC---AAATCGTTCTCACTTGCC-CAGTTACACT--CC-G-T--ACTAGCTACGAGC-CA-AGAAC-TAAAACTCGGCAACCGACTC-CTATAACTC--CCTCTTCTTCTAGAGAACTCAGGTA-CGCATGGCTAAA-AAACAATGCAATATCATTAAAAATACTCTTTAAGCCT-G--ACAAAGAG-AGAC-ACAATC-CAGGTGACTC-GTTACTCCATGTGAGAAC-TATTAAC-ACTTACCA-AGGCTTAAAGTAC-TATCC-ACCTATTCCTCAGGAGCCGAC-CTAAGGTCCACCTGTATTAAATCACCCATTTGGGAACCCACCGGGAACTCGCAGTTTCAACATACCGCAAACAAAGATTCATA-AAACAACCAATAACC-TCCAGGAAT-C-AATC-CA-ATCACTGTACAAACAGTCCACACGA-GGCGAGCAACTATTGATGGTATGACACTACATTTCGAGCATCGGGGCAACAGGCGAAAGTCTCAT-C-ATAAGAACACT-TACTTTGCC--A--CCACACATTGGAAGTAATTCTAAGACTCCTAGTAAAGGTAGGATC-TGTTGTCAAC-CATTGACGTGACGTGAGTCAAAGTTCAACAACACCCACCCGTAACACAATCTCCTGTAAGCCTACCTGTCAGC-TG-G--CAATCGATTC-A-TC-CCTAAATAGAAACGCTAGGAACTTCTC-ATCAGCTATAACACTTTCGAATGCACA-CCGGGTGACCTATATTGGTTAAGTACTCATTTCGCCACAAC-GCCTATAATCCGTAGTAATCTACGTAACTAGTCATCCATATACCGATACTGTTGTAGACA-CCCGACTCTCCCTCTCTCAG-CCCT--CCAGTC-ATGAAACTAAACTTTACGAGAAGCCGACTCAATTGTTCTTCT-CCTACTCATGTATAATATTTCTAGATCACTCAGTCCCTTACAAGGTCCCCTTAAG-AACAGAGCCATTCACACCTCAGCGCAAATCTTATGAGCCGATGGCATTTCACGTGCAGAGAACCCGTTACTGCAAATCCCGCCCAAC-TCCGCGGTAAGCATAAACCTCCAATCGTCCTCTCACAGCCCCTCCAAATCACTCCAAGCACGTTAAGCCAGTAACGACTTGTAACACTTCATCTCTTACATCGAACCATAAGAAATTCACATACCCTTGTCAATTACCGGAAAGACACGCACTTGATAGTTATATAGACTCTTAGACGAATACATTC-ATCACACATTT-TTATTCCCTCGCAGATCAACCACTCACCCACACTCCCG--ATTAAAAATTCCATACCATCC-AGAGACCGAT-CAACCTG-ACGTGATCCAATCAA-C-CA-CACCCGTCC-CT-GACATCCCTATGTAGC-ACCGACCGGTCAACA--ATA-GGATGTTGCTTATGGTTAAGGACACACTTGAATAT-CTGCCTCA-CAGAGGTATACTAC--TTGTC-CATCGGATAATAGATGTATACGAATAAG-AACTATCACCACTG--C-GACCAACACACAACACATCGATAGAAAATAAACACCAGACATACAAGC-GGCGACAAATTAAG-AATAGCCGACCCCGATAAAAACG-GCGAAACCAA-AATGCTTCTTCTCTGATTTATTTCGGCATCAC-ATGACC-CCTAAGCTCTAAACTCGCTATTAAAGTACCTATGAACTCCGCGAAAATGCTCCCGCCCCGCCCAAACAAC-CTCTTTCC-T-GCAC---CGGAG-CC---CCCCAATCTAGAAGTTACCGCTGGGTAGCGACTTTAGGGAACCCACCCACTTACCACACACATGCC-ATTCGTTGCCCACCGTCGTGTTTCCCAACCACCAATCT-TCGTAAAGTGGAATATCGAAAGCGGCACTCTATTATTCCTCTGTT--TCCGCATAAGCAG-GA-GGAAACCACGATAAACATAATAAATT-TGTTGCACACCACCCACACAATACCCTTCTCACCCGACATTATATGCCAGTAACACC--GCCCTAACAGACCCGGGCTAGAAAATACTTTC-GTAGAAGCA-GACGGCAGGACCTACACGGC-C--ATCGCACCCTGA-AACCCCCCAGC--TGCTTTAGTACCTAT-ATTCTTTGCGTCA-ACAACCA-AACTATACGT-ACCTTCATTGCAGAAT-TTCC-G-ACAGTGCAAC-CTGAGCAGTAAGAATCTTGTAACTCTACTGAAAACCCAAACGGCCATA-AG---CC-GCTTC----CCACAAGATTTTTTACTTCCAAAAT--TC-GTCATC-CCCTTACAAATTCCGCCGCA-CCTAACAC-AC-ATAATACAT-CAAATCGAACCGTTAGCGACAAAAGTTG-ATCGACAGTAGAAATAAA-AAG-TACT-TGCTCAACATAGATTC-C-AC---CCAGATAGCATCACAATCCAGTTCTG-AAAAACC-TTCACACCTA-C-AGCTAAAAACTGGTGGTAGCTGTTGGATAACACCAAAAAATTTGACAG-GT--TTACCCGCGTAATACAGAC-AAGAC-ATTTGGTAAACAGACTGACTAATGCTTTTC-CAAACA-TAACA--AACTAAGCTTGATGTAACATCTTTCT-CAACCTATCCTTTCTAGCTTCTTCCAAGTAACAAATTTTGACCGATGCTATCTCCTCTTAAC-CAGGACCGAAAC-TCTTAGGA-ATATAAGCTGAGGAAAATA-ACAATC-TTACAACGCATCTGCACCAAGGACCGAGCCCAGCGAGAAAATCGTCTTCAACACGCACCCCAATATCATAGACTTTATAGCAAGAGACCCTGCGCCCTCCCTTATTTATAACTTTAC-GATATCATAAAATCAAG-TGTTATAAAACAGATACC-ATTTTTACCTCAAT-G-TTTCAGACTCA-CTA-CACCACTAACGTTGCA-ATAGAAGATCTTGAAAAAACCTTGCTTCAACTGGATTATTATTCAAAATC-A--ACATA-A-GCCCA--TTGG-AAT-CG--CCTACGAAGACATATCT-C-TATCCA--C-CTCCGCCGAATGGAAGCCAACCTAAGTAT-CGCCACC-TCATACCCCTAAAAATGTTAAACATTCCACCGAGGTTGCCC-GCTCCATTCTAAGACTTTTAAGTT--TAGGCCCATG-GCAAGGTA-TGCTCGCCGCGCCAC--A-ACGAACTTGGACTATAACGGTCACCGACG-TATATCAAGTCAGACCCA-ACG----ACTTCGC-GC-CCCGTATTGACAACGGAGGCAAATTCCACGA-AGCTCAGCTGAGCACAC-AACAGTAAACGCATA-CGCTTTATTTCAC-CTACAAATTCCTGGT-CAC-TC-AGG--ATATCTAACCACATTAATGCCAAAAATACAACCAT-ACTCC---GACA-CGTGAAACTATCCACAACCCCCTCCTTTCACTAACCG-CCCTTTAATACCGCTTTAAAGGCTACAATTAGACAATCTACAGCATCATTCAAATGAGCATTCTC-CCTGTAATACTTTCGA-GCTCAACTACCC-----CGTC--G--ATCCGG--AGTGGC-AAATC--CTTCC-TTGCA-CTA-ATACCATTCAAATCAGATAACA-CAACTCCTATAGAACAACACTA-AG-ACCTTTAACTTAAGCACTCATAACCTTAACCTTCACTGCTCCAAACTTCAAATCTCTTCCAGGACAGGAACTATTGTAAACGAGAGAGTTTAAAAATC-CTATTATCTTAAAGCACGACATTATCCTGCGCACTGGAAGTATCCAACTTGGATCC-CAA-AT-TAAG-A-TTCCATTCA-TGA-AGAATACTATAA-ACCAACC-ACG-AAGAACAATA-TAGGCA-T-CTAAAAC--T-CCCGGTCAACATAGCC--GAAACGTCAAGGCCTAACAGACCTTCTCCGG-AATTCTCTACTATACCAAAATCCCCTGATCTAAAAAATTAAACTGAGCTTAACA-AAATGGCAATGT-CTTAGACTTACC-AAC-ATCATACAACTTGACCAGAATACTCTTTTCCACAAAACCGATGATAATCATCGATCAGTCCACCCGTTAGACAAACTG-AAACTCAATAAG-A--AGCAGACC-GC--ACCCTACTAATGGCTAAATTATATAC-TAAT--ACTTGATAA-CTCA-TCTTTTAACATCGAATCGAC-CCGTTGCAGC-CCCA-AGAATATCACTATCGCCAGTAGCTTCGCCAATTTCTTAACGTTG-TAAAGT-CCGAAATTAGTATTTAATGG-CATC-AAGGTTCGAAAC-ACACTCGCTATCCGC-GACAAAGTAAAC-TACTCACTAAATCCT-CCAATTATGTCTTTCGACAAACCAGAATTTACG-GTGAATCACT-CAAGCCCTCCTCACCAACGC-TAGATAAAATCGGC-ATCTAGACGGA-GT-CTCCGTTTAATAAAATAACGGAATTCGACGTTCC-CCAC-ACATACCAA-GGGCAACGCTTTACTCCCAAACATTGAAGTAAAAATGGGGCTTCTTAGCATC-CATAGCCAGGGCGAGTCCCCTAATACTCCTAGTGAGAAAGCGAAAGGTAAAGCGCGTCAATCTAGCGAGAATACCCAAACTGCTTAGTCACGATCTTTGTCCTA-CATAATTGACTCCTATGGCATCTTAGGAAAGCATCCGACAGAGGATCTATCCTGCTGATGCATTGCGTTCGTATTAACCAACAATTTACCTGTGAGC-TCTAATCT-TATCATATGCC-CC-AACA-AAACTTAATC-GGTTGCTACATCCA-GCAATAATTTTCACCTCACATTCA-ACTGTCCACGTATGCATTAGCTCTATCCCATTAC-CCT-AAAACCTATCGCCAAT-CTT-AAG-CGTCAG-ACTCACGTAACCTCAAAACTGACTACCTCATTCCAAATCCAAGACTACGATCCTCCCACACATGTTAG-CCACGACGATAATAAGCGAACCAACCGTGGAAAATCACATTCGACTTTGTGC-TC-ACCTCTACTG-ACCAT-CGACCACGTACACTTCTGGTCGTCAACTATGC-GGGCCCTCCGCCAGACTTACAAG-CAACTT----CTTCGCAGAAACTGACATTAAGAAAAATTAAATACAGGCGCCATACGCCCCAAGTTATCCATTCTTAGATTCACTTTATT-CCATCAAACTATGTATTCGCAACTCACAAGGCTC-CCAATGCTGAC----TGACTACTACCCAGACCTTCATAT-AAACTA--CGTGCAA-CTG-TAAAATGCTTCACGACCATTCGCTTTCCCCACAGTCCAGAACCTGA-CC-CTA---AGAAGGTAATATAACGCCAAACCAT-CGCTTCAACAAAAACATAATAG-CAATTCCAAAATGA-ACA-GGATAGATAAGTTAAAGGTTTACAGTCTTCCCC-ATAACTATGTTC-CTGATCC-TTAAATATAATGGGTTA-GTA-C--AAC-GCAAACGCT--CAAATCGGATCATTCTTGGA-TAATA-ACAGTCCC-ATATAC-T-TTTCA-GAAATACACACGTGAGGGAAAAC-TA-GAATGCCGT-AATATAGACACTTTC-T-GTCGAG-AAGAGGTG-AT-TTCG-TTCATCCA-TCCTTAAAGTTACCGCAAC-CGCTTCTCTAGCC-GATTC-ACA-G-AC-C-CTTCCCCCCAAGCTTACT-AC--TCC--CTGATCACTACCATTAGCTAG-TGTACTGTCAG-GCCGGATTCATATAGTCACCACCACGCTATGGTCA-CCCCTAACACTAACCTCCTTAAAATACCAGCAATAGGCCATACTCCACGGATAGATTACCTAAAAATCAC--ACCTAACGCGAACAAACTAGC-CTATTATTCTC-ACAATTGATTATTACACCAGACATCCTAAATTTCTTCTTGCTAATCCCTATCTCTCGACGCTAGCAAAGAGACTACTCC-TAACGTTGTTTCTAAACACCTCGCCAATCTCT-TAAA-CCGAAGTTTGACAATACCGTTTCTCAGAACCTCTTCTAAATATGGCAAGCTAAC-TTTCATCTTAGAGTTATCCAGG-TAGTGAAAAGTATACGTTACGTGCTATATCAGTATCCGCAACTATATACTCAAGGATAATCTTCC--AGCACACACTTA-TA-TA-CACCTTGATTTATG-TA-TCATAAGAAATACGCCACAGTCTTGAA--C-TGAATCGTGTGAAAATTCTACT-CCATGAATCCAGTGACCTTACGAATATCGCGAAA-TCATCA-CTGAGCGTTCT-TAGACCCTTCTGGCACCTTTATGTTTTTAAAATACT-AAACATGCAA-ACA-ACTTC-TTATTGTAAAGCCTAGCACTTCAGC-CACGGATCCTCCG-A-AAACACC-AACT-AGATC-TGA-G-T-GACAACATGAACTCC-CCTTA-ATAGCGATTGAGCACTCTCATT-TA-CC-G-A---CG-CAACTTTTGGATCAAT-TTATCC-CCCCCT-GATTCAATGACC-TCCGCATCTAGCCACCGCT-T-CTA-A-CA-TGTGTTAGTGTACAATAAA-AGACCCGTACAACCCAAAAAATCTGCTCA-CTCCGCTATCAGCCCAAGGACATCAGCGTCATGG-AACCCCATA-TCCAA-ACCTCAA-TATAACACA-CCCCTCGTCTCCGTAATAATTATAG-TTACTCGGATTGTAAACATCCCTCGCGTATTCTTCTCAGTACAGTAATAGCTTATCTTCCCGATGTGCACCCTCCTGCGCGGTTTACCTA-AACCATACTCATCTTCCACCAATCTGGTACATCAATTCCATAACACCCGATG--AAGTTTTAACCGCCTGAGTGTTATTCAGCTTCAGCATGTCAAACT-TACCAAATGCC-GATTGTG-GAGCTCACACATTATTGAACACAACTC--T-TACTG-ACTTAATC-A-TTCTACCCTACCGCCATAAAGCCT-GTTA--CGAATAGCGCATTTCTAACATAGAATCACCAACCATCTAACGATGATTCTCACCCCATACCCCATCATCTAAGTGCCACCCTCTAGCTAGAACAGTTCTTCCACCGCAACGTTCCGGGTCGTACACTTAACGTTACAAACGAGCCTAAA-ATTC-AGGGTATACTCTGCCCTTT-GAATAA-CTAGTCATA-TTCAGCTCTATAAACCAAAACATTATTTTAACTTTCCCTAAATTTACGTCAAATAAACATT-CCTTAT-CGCTTACCTCTTTAGTTACTGCTCGACGT--A-CGCATGCGACT-ACA-CTACGGGGTACCCAGAAAAAAGCTCAACCACGCAATAACA-AACTGTGCATGAACCGTTAAACTTACTAGAAATAGGTCTCATCCCTGCCACTAGGATCCAGA-CAAT--CAGA-AT-AGGACACT-T-CCCCTAAATCAATT-GTAATCTCCA-G-AGGACTAGC-CAAAT-C-ATGCCTCGTAA-AACTAAATCATAGGCAC-CTCCCTCGCACC-GATACTCTCCCTAAAAAACAC-TATC-CAGATCCCCAAC-TACATAA-AACTAAAATCAGAC-GCCTAAAACTAATATCTTCTCA-ACTATAAAACCGATGAACACAC-ACCAAT--CACTA-C-TTT-CCTTG-A-CATAGCCCTC--ACCAGGAATTTTTAC-GCAAAT-GAA-TATCT-CTG-TACC-ACTGACTTAGAACGAT-TCGA-C-ACAAGAACACAAACAAATTC-ATT-TAAGTTTTTTCGCTTTAC-TTTCGTATACTGTGAAGAC-TCCAGCCCTGTAGACCTTTCCCAACCATTAAATCTAGTC-TTACTAG--A-CAGTATAAGACTATGGCA-TACGG-TCCTATAAATACTAATTTC-TGAACTCTATCAA-TTATACAACTC---ATATT-ATGAATGT-CCTGAC-CAGTAC-CGA-GCCAGGCG-AGCGATGATCAGAGCGTTAAATATCTCAAAGCCATCCCAATTAACGCCCAAGAACCCTCACTAGAAATAGTGCCCCGAG-TTCTCCCCACTATA-AC-AT--CACCAC-T-ATA-CA-CCTAGGGATCAT-ACATATAGACAG-CGCTCAATGGTACTACTAAAAGTCGAAACGTTT-TGAGAGA-ACAAA--ACCA--CATGTAATATCCCGTGGATTC-TAAA-TTTCA-GGCCAACAG-TGT-ACACAA-CACGGAAACTACCACTAATATCTTAATAAAAACCCCAGAATTCAGGTCTAAGAAAACATAGATACCTCATCCTCTCCCTCCGTGCT-TAAG-GCAATTA-GCATCCT-CTGACATATAGAATTTAATATTACT-AAGAAAACGACAACCCAGGCGATTCTATTGTACACTCGACTATACAAAGTGTTTAGC-TGACTCTCGCGCACCCGTATCATTCTCTCATTCCCCTTTCTCAAAAACACGTGCCCTTAACCAGCCATACGACAAATTTTCCAACATTTAATCACAAAATGATCTAAAGACAACCATGCATT-CGAAGGATAATCCCGCGCCCTCATTAACCT-TAACGCTCTCTTACTGATCCCACGACA-ACGACACTCTTA-C---TAAACTCAACTACAG-CCTCTACGCTCCAAGTTAATCCTTACAAGT-T--TTACAAAGTAC-CGAAGG--TACA-ATGACCCAGC-CTTTGGTATTCCACTC-AGACTCCTCCAAATACATATT-TGGA--ACAAGCATCAGAA-ACCATTCCCTA-AATCTTTCTACAGGGAAAGAACCCTAACAGG-GAATCAAACACC-TCGAACTTTGTCCTCT-CTACCGACGGCACATCGAGCATCCCCGATACCATACAGCACT-TGCCGCTCATACCAAAGAATCGCCGA-ACAATCGGGTCCCTAATTACT--ATTTAG-ATCCTCCAACTCTCTCTAACTC--C-CAATAGTAGCC-TTAGTGTAACGACTCAAATAACCATTATATCAAT-A--AAACAGTGCCATAGAAAAC-TATTATACCGCGTTTCAAAACCTAAAACCCTGGA-CATCCGCGCCC-ACAAGAAATTCCAGTCC-A-TCCACCTAACAAGTCTA-CAGA-TAATAGCTAC-ACATAACAGGCACGACATCAAAAACAAC-CATTATC-CGTGCATCTC-AGATA-TTCCAGCTATATCAAAAAATGAATTCCCGTTCGCACTAACTAAAACCACCTGTTTAATTATGTACGGGCTCCATTCACTTAACAACAGAGAC-AGAACTCT-CCAGAAGGTCCT-TGCTGGAATAACATTTCAAGTCCC-TAACC--ATCCCAACGAAGATTCATCAGCAC-TCTGG-CTGAGAATCTTCGAGACAGCATT-TAATCGGCCGTCATTGC-ACATGCGTTCTATCCATTTGACACAAT-CACACAACGACAGCCCT-G-AACTTCCAATTCCATAAACCCTCGTGTGC-GT-TAC--CTTGTGATACGAGAAAAGTGAAGTCATCG-TCA--AG--CCGCCCT-CGATCCGCCAACATCGGT-ACCGAAG-C-ACACCCACACCTCACGGGGTC----GCG-A-T-TCTAACCCATCTCATCCAGCTA-ACTCTAAAAAGTTTAACAAA-TT-ATAT-AG-ACTTCC-CGCTAGTAACAAGAAAC-CGTTGGAA-T--TCATAGATTTAATATCCATGC-CGA-TTGTAAATGACT-ACTC-ATAACCTCAAACAAACGTTTTTTAGCCCTCGCTAGACC-GTTAAAACG-CCTATGTCTACTA-AATCCAAACATTAA-CTTAGACTACGCTAAAATTGCATTCCCACACAGTCACCTCCTGA-G-A--AGGATTGCTAAGTTACGG-AGAGTACG-T-GTCCCGTGT-AT-GAC-TCCGCAACAGAGC-C-TA-GTCACCC-TG---CTC-CCAAATCACCCAAAATGAGACA-AAGCATAATACAAATTGCTTTCAGATT-GTCAATACCGC-ACCGAAGAATGCTACTCATAACCTATCTAAACTATTGAAACTTGGTGTTTA-TGAAATATGGCTCAGCGACCCCTAATGCTAAGGCCAAAATGTCT-AAGACTTGTCAATCAAACGGCATGAAACCAAT-CCC-CT--CT--TCAGA-ATGA-CC----CAGTAGTTTAGG---A-CAG-ATCC-TTGCC-CTTTC--ACA-TCATACCCTAC-CAACTAGAAACAGGAGGCAGAGTTGTATAGCAAACACATAAATTGTCACCCAACCATCCTTCTAACC-CACTAGAGGATACT-GTG-CGGATTCGCAAATCCTCCTCCCCCGGTCTGTTAAACCCCTAACTT-AGAA-CTCTTCCATTACTCA-GCTCCTCAAATACGACC-ACTAAACCCCTTATCACGTCTGC-AACTTACCGCAATGCC-GAC-C--ATCCGGCATATCAAATTCAATAA-AG-ACTTG-CCCACACCAAACGCGTTAACCTAACACCCTCAACACCTTGCTAATTTCATTAACGACACAACAGAAGCCTCC-ATAGC-CG--CT-AA-CGA-CCCCCGGGTAAAGTA--GCTA-T---CGA-CG-TG-G-GAATAAACACTTGTACCTCAAAAAAGTCCAAATGCCTAACCATTTTCAG-TACTTATACTGACAACAAATAGTACTCCCCCAACAC-CCCTTAAACTAG-TTACTAACTTCACGCATTATAATCAGAAATTCTGGTCAACAGAT-CTCCATATCCTAATCCGTATGTAAAATGCATAC-CGAACCATGTTCTTCAGGGCACAGTAACA-TCAGAATAACCCTATATAACTAACTTTCCGTAAGTATGTGTGGTCATCTAAAACCATAACTTTCAATGCGTGTTACTAAAACACTTTAAATACGCA-ACGGATCCAGGTTACAAAAACACGGAAAGAACGATTCACTACAC-AATCACGACC---CCTC-CATGATCAT-TTTCCCACCCGCACATTCAC-G-C-ACCCATATT-GCAA--TAGTACC-G--TTA-TCAATGAT-TTCTCAAAACGACCTTCTTGATAGTCAAATTACAGGTCAACCTATTAAATAAACA-AATACAA-TCCAAACACGCTTCACTCGCTT-CAACCCAGCAAAATTACGAACCCCGGTCCAAGAATCTCATTATAAATTTCCACCTCTCCATGTTCAACTATTCTAAGT-TTCAATGAT-GCCCTACACCT-AACTTCCAAGACTA-G-AAACCCCAAACGAACC---C-CCTGATCCCATTTACCTCAGAATACAATTCGTACTCCGAAG-C-CGACAAAGACA-GCCGT-ATAGC--CCAGCTCGCAGAAATGACATCGGCGAGGGTAAATTGACAAGGGTTC-TTA-CATTTCAACTTGACCCC-ACACCCGCCAATTC-CTCGC-AATTCGCGCACTCATATTCGTACACGTA-CCAACAGGGTTATATCCATTTGTTGAAAAC--CAGCTAC-TAAAAATAAAC-CCCTA-CAC-GCAGAGACTTC-TATCGCC--CACCACCTACATAACTGGAATAACAGGC-AAATGCAAAACTCTCC-GACTTTCGGCATTCGC-TGAC-C-CCCCAAATTACCCATAGAATTCTCACTCTGACGAACACGCCGCGCCACGTCC-GCACATAACTTAGACGCAATAGGCAGACTCTAG-CTCTTACATACAATAACCATTGTAATTAGATATGGCTCAACTG-TTCTCGAAGGACA-AATATCGCCGTTACCCTATCGCA-TA-CTCTACTC-GATACCACTGGATTAT-GA-AAA-GTGCGCAACGCAATATGAATTAAACCGACGTATTCTCAAACTCACCCAGCC-CA-CAC-TAATAGACAAGC-GCAAATT-T-AACCTCTTATTCTC-CCTTAAAGCTGGTTAATCCCAGCAATAGACACTCTAGTAAACTGAGAACGACCAAGACCCGTTCTGGAACAGAGATCCACCTTTCATACCGCCATAT-TCAACGCCCGATCACAATC-GAGTAGTTTGCACTCGACGAGGTAGACTATCTTC-CA---CCTCCTTGTTTAAC-ACCTATTATCGATGCCGTACCTAAAGGTTCCATAAACTTCAGAAATCCAACAATGAC-TCTGTATTCCAAGCGCTGAAAACAAGATAGAGAATCATAC-CCAGTACATCACTCTAATACAGAGCATACGTCGAAC-CTAT-AGCCAGAACTTCGGAAACACAC-TCCTTCGAACCGTTTACCAACGGACGCGTTAACT-AACTTCTT-TCTAAAACCTAC-CAGCCTCTCCTATCCACTAATA-ACGATGA-GTAT-AGATAACTAAC-GCGAAGATGGCAAACACACAG-CCT-GCGGCATTTAACTTGGCA-CACGACACCT-GCAATCTATATATC-AAGCTCCAACGTCAGTCACTGTC-TCCGCTCTCAGTACTAATCC-CACACA-AACT-ATCAC-TTACGAAATAATCTT-TGTCAAATACGTTGTACTGACTGGTCTACCTTAAAGTCTCTAACCTAAGTTTATCCCATCCTTGGCAAACTCGTTTTATACAATTAACCCCCACGAAGATGTGATTACGTGGCTTCGTATAATCTCACCGCAAATTTTCAA-AAAAATCCTGTCCTAGGAAACTTATA-TCGTGAACATGCCCGAACTC-TCTTCGAGAACCCCTATTACCAAAATTAAACACTAACAAAAGTTCATATGATTGTATTCAGACATCCCTTCTAAACCAAGTAAAAAACCCACTTCCATCAGGGAAACACATTTAAACTGTGACCT-ACAAACCTTCATACTTGTCGCTATACAAGAGACTAGCCATACACCAAACATCGGTCGAC-AC-CAACCAATCCCTACATGCCTTTTAACCCACCTAAACTGACGACAACCGCACACGTCGTCGACACTTTAA-CTTGTCAATCGCAATAAATACCGAGCTGC-ACA-CCGA-TATCA-A-A-GAAAC-AGAAGAG-CGACTCGG-CG-TCGCT-GTCCCCCTAAGATT-A-AA-CCCTTTCTTCTTTAGCTTAGCTTCTACGACCCCATAATGAATTGTTAACCTACTTTTCATTACAGTTAACCGAC-AAATACTAAAGCGTAGAC-G-ACCCGCGCAAGTCACTCCCT-ACACA-ACTCTAATACATAACCT-TTTAAGCTACCGTCTCACAAGTTGCAG-TACCACCTGGATATCTAGCCGGTCTCCTC-CGCACACAGCTCGCTCTAACTCATAGTGAATCAAAG-TACCGTTG-CCT--CG-CCGAACTC-CACAACGGACGTTAAGTTTAGCTC-C-CC-ACAC-GCTGACCCCACCCG-CATGGCGCCAAGCTTCAGCAAAAACGC--ATT-CC-AC-CCC-CATGT---TCAAAATAGT-TAACCA--TTCCT-AAGCATAAACATTCTAATGTGCCAAATGCATATTTCAGAAAAC-GAC-CCCCAGAATTAGATACAACAACATCACCTTTGAGCTTCCCATAAAACTAGTCTCGTCACTTAAGGAATAGGATCTTTAACAGATCTCCCACGTCCCCAACCACCCAATAATCTGCTCC-AAAGAGCACCCGACC-T--ACCTTA-AAAATTCTGGTGACATTATAC-GTCTAAACGAATGGCCG-GACACACG-GTGACAGTATGCAGAATTAAATCCACACACTACC-----CTTCTCAAGCGTA-TAATAAAGTAC---CG-TTTACACCCCCTTCCACAGAGATGGAGCTCTCCGAGCCATTTATAATTGCAGTAGCTTCATCTTTCCTAATGCTACACAGTATACGTCGACGA-TTACACCCAATAAGTCAATGGAACCTGCCCTCAACGCCTGCCTAT-ATTCACTCCCTCTAGTCCAAGATGCACACCTAAAAGAC-C-ACCCCAATACAGTCACACTACTTAATACCACATACCCACC-CATGATGTACGTTGCCC-CAATTACCACACTCTTGACAAGTGGTGACAAGGAAGTATGTAGGCCATTAAATTTAAACTATGGTCACTCTACAGTACAATAGCATCCTGTTACC-CGAGTGGCATATTACTAATCGTGTAGCCTACCCTGCCACTCTCGATCTGTATCT-AACAAACATCCATATCACTTCCGTATACACGGTATCGTCCGCGCCAAACTCCAACTCATAAATTATGCCTAGAACACAATATTACTCTCTTTTTGGCAAATC-GCCTCTTAATAT-C-A-CCATCT-CATCTTATAAGTCATTCAT-CAGA-T--TATAGCTACGTA-CTC-CTC-GTATATCTTTCACTC-CGCGACCCCAATG-TAGTATACT-AGCCTACAACAA-ACT--ACCT-CTC-T--ATCATTATCGCGAGGAATACTGTG-ACAGTTATAA-AAAC-AC-ATACTGA-AACTTG-GTCGT-ACCTTAC--CTTG-GGCC-TCAG--A-CAAACAACAGGAACTGAAGTTGATCCACAGATTACCGGTCTCAATGTAGAAAGC-CTACTAC--TGCTCAAATACCTTCAGTCGAAACAAA'
    ' ||  | |::|:||||  :||| | ||| |  || | |  | |||||:|::| ||  ||:|:||  |: :|: ||  |:::|| ||| | |||:||:: || :|| : | |  |:::|   ||  | : : :| | : ||:|:| |||| |:|| |||  ||:|||||  :| ||| :|:||| :| | : :: | | ||| || | |||| | || |  :::|  :||||:| || | :|||: |||: || ||   ||::   :  ||: ||||  || |:   | |::|::  : :| : |||| :| | |:||:|  |::|||||:|||:|:| |:| ||| | ::| : :|| :|::| | ||:||  || ||| ||   | | |||  ||   |:| :| ||||::: | | |  ||||:|| :| ::||||  | || ::  ::||::  :| |||   ||:::| :||:::  ||| |:|| : | |  ||: |::|||||| |||::|| |:|| :|:|:|: || | |  ||| || |||: | | |: |:|| |||:|: |: ||| | ||||: : ||: |:| || |:  |: |:| || |||   | |:  ::||  | || || || :|| ||::::::|  |::: ||:| ||  |||:|||| |: | ::   | |||| |:: | |: ||| :|: :|||   |   | | || :| :||   | || | ||||  ||| |  ||: |::||  || | |: ::: | |:||:  :|:|||: |:::|:| :|:| |:||:::|:| :||::| || ::||: |:| |  ||:|  :||| || |:::| ||||||| || |: |  || ||:|||:| : |:||||:::|||   || | :|| || :|:::|:||: :||| | ||::||:|    |||  ||| :::|:||:||| :||| :|||:|| | :|| |: |:| |:|:::|| |  :|::|:|  ||| ::||:| ||  | |||||  ||: ||||: ||    ||  :  :|  ||||| | ||: | ||||  || |: | | |:|| || |:|||| | ||  ||:|| |:|| ||:||| ||   ||| : | ||  |:|: | ||||| |  | |  | ||| | | :||| | || |  : || | || | | |::||::|| | |:| |||| :|:|| :| | :||| |||: |:||::|| :| | ||| | || | :| : ||| | | |:| | | || ||: |: |  ||| :  | | |  ||||| | |||:  || ||   |:|  || || :| | || |: || | :::|:| | || :|: || ||| : |  || :|: ||||    :| |||    : || |: :| |::||||  | |:| |  ||| ||||| |  |:|  | |:: ||||| |:: || || :  | || | :|| ||||::|:| :||:  | :|  ||:||| :  ::|||:| ::|:|| | || :| |:|: | ||| || |   | |::|  |||||::  || |: | :|| ::||  ||| :|:||||  :|:| | || ||| | |::  || ||| ||| : | ||  :|:  |:||  :||:: || | | | :::||:   :|:| :| : | :| |||||::|  |:| :||  | || : :|: :   ||  ||||||: ||| :  ||| :|:||  || || | | :|| |:||  | :||| |||| |::||| :|::||| ::|:| | |: :|| || || :|::|  | :|:|||   |  :|| || ||| |:|  ||:|:  :|  | ||   | | || |:| | | ::||||   || ||||  || | |||||:||:: | ::|| || ||  || :|:||| | ||: | | | |:|||| | | || | : ||  ||  ||:|| :| |: |||::: |:: : ||  ||: |:|| | ||::||  |:| :|| |:| |   :|| | :|   || | : :||||   |  | |:|| |  :::||| | ::  ||:|| |:| |||| :|  :| | || |||: | ||:||:|||: | :|| | |:| :||   ||: || |  ::|:| ||  |:|:|:||:|||: : :  :| ||| |||: ||:| | || :|:|||:|| |::|  :| || |: | ||||  |:|| : |||||| | || |:||: |||| || |: |||||||:| |:||:  | :|| :| |||:||| ||||  || |:  :: ||||||  || | :|  |:::| || | || |:||    : | ||| | | ||| | :|:||    |||  |:||||:|  |||||:| |   :|| |:|: |  ::|:||  || |:|: :| || | | :|  ::|: ||||  | ||||||||:  : ||:||: :| |  |:|   :|| ||| |||::| | :||||:|| || |:||:  |:||:|| || :||:|:||:| || ||| | | :||  |  || | || :|||:|| ||: ||  : | |: |:: | :|: |: ||| || |:| |||| | :||  ||||:||:: :| | :::|||| | :| | | | || | :|| :||:|: |  |:| |:::||||   |: ||||  |:||:| ||  :|| ||:| || |:||: ||||  |:  |  |||: |:|  || | :|||: :||| |:|| ::| |: |:| ||| ||  | :||  || ||: ||||:|:||  : || | | || :||||  |||  |||| |:: |::| |:: | || |||::: ||  || :||   |||: :|:|| || :||| ||  |||||:  || | :  |||  | |  : | || | ||| ||||||| | ||: ||: :|  :||| |||  ||:|:: : | |  :|| : ||:|: || ||  | |:: |::||::| :|  |||::   |: | :|:  |: ||| :  | |||:|: || | ||:  | ||   :|| | |:|:||||| ||||  || || |::| ||  || | :|||  ||  ||| ||:|| ||| ||:| : |::|   ||  | |  ::||:|||||| |  : || | :|:|  |::  ||| ||| ||| : :: |||:| | |  :||:||||: | ||: | |::|:| |  |  || |||: ||:| | |||| |: :||: |:|:|:|||:|   || ||:  | :  |: ||||| | ||: | || | | |  |||:| ||   :|||:|| || :|: || |:| | |||| : |||   :|| | ||  |  |||| |: |||: :: | ||| :| | ||: :: ::: :|||| | ||   | ||  |   || :|| || || || :||   |: |: : | |  ||:|||:|||||: || :: ||| | || ||: | || |  |||: |:|:| | || ||| | :: | | | || |||| :|:| : || ||:::| :||| | || ||:||  |:| |:: |:: :|| :: :| |:|| : | |:||| ::| :|| : ||:|| || | | || | :|:|: |: :|||:||::|||   | | |: |::|| |||||| |: |||| | |||   ||  ||:||| | |:||||:::||||:| |::|| | :| |||| | | || |  |:| ||:| ::|||: ||:||  | |||| :||   ||||    |||||::|:| |||| | || ||||||  :|| |||  |||| : || ||:|| | | : || ||:|  |  :|||   | || ::| |||    | | |||:  ||:: | | :|||:  | :|| ||:| | || |  | ||||  ||  |||:: : |: | |||   || ||| ::||:|| : |||| |: :||:  :| || |||:|: ||||| |  ||:| |||  :| ||:::|:|| | :| :   :|  : |:|:||:|:|: | |::|: || | ||| ||||||| :| |  ||||:||  : |:  ||||| | :|: |   |  ||| ::    ||:||| | || |  :| ||:|  |: |:| | | :| || :|||| |  | |: |||: | |: :|||:|| |  ||| | ||| |||:|||| ||  :|| :||:   | ||||:| ||  | |:||||: |:||| || || || |||  |    |:||  || || |  |: :|::||:| || |: :| :||| ::| |||||||| |  |  | |: |||  || ||: |   : |:| : :||:||||||| :|: :| |||:| ::: | || |:||  | :| |  |||  ::|   |: |:| || :|: ||: |||:|  | | : ||||| |||:  |::| |  ||||| | ::|||:|:  | | :|  ||| | :|:| ||  |::: :||: |   || :| ||:| |:  ||| |||| | ||||   ||:|| ||    | | ||||  |:|  | | |||: :|:|:| : :  : :|||| |||| :|||| | |: |::|  ||||  : |||| | :|:: :|:|  || |:  ||:|| :|:  | |   |||||| :    :||| |||:||:|: |||:||:  :| | :|:||||: :| :|||:| ||:|||:||::    : |  ::|||| ||  | ||: ||  || ||: || :  |: ||:|  | |||  || ||||  || | :||  :||||:: :  |:::| ||:|  | :|||: | |:|  ||  | |:| ||:|| |  ||   : ||||: ::| | | :||  ||||||| |:| ||:  :||: | || || ||:: || :||:::|| | || |:|| :| |:||  | || | ||||| | :|  |:|||| :|:|| :|  |  :| || | :||| ||| ::| || :|| ||   || |||||     ||||:|:      :| :||:: ||  || ::|  | |||::|||: ||| | |  | |||:: || || ||  || :| |::|: : :||:    | |||:||:||| |||  | | :   | ||| | | | || |:|: |:|| |:||:| | |    ||:::|||| ||:|:|  |: |||:: |||:||| :| || ||:| | |||  :  |||  |||::||  :::|   ||: |:::|||| ||| || :|  :|: ||: :|: ||||: | :: :| |||||:   :|| :| ::||||:|||:| | |:::|| |:|||   :|  |||||  ||::||  |||:|  ||:| |:| |:   ::| | |||  |:||: |::| |  :| || |  |:|||| || :::| |:  |||  |:|  | :  || |::|::| |: | :  ||| || ||| :|||:| :|::|: |:||::::||| | |||  | |: |||:  ||  ||:| |:||  | |::|||:|  :|  |||: | | |  |||| |:||| | |: :||:::  ||   | ||  || ||:|:| |:| :| : || :||| ||:|| | : |  || |::| :  |||: ||||: |:| |||||||  |||:||| | |||| : || |:::||   |||||:|| |::|:  :: | ||:| | | |  :|:|: :  ||||  ||:| ||| ||  ||  |||::||| |:|| | | :|||     || | | |  ||||:|||:|| |: |||  |  ||| ||  ||||| |::|| :|:|::|:||| : |:| | : ||| ||  |:: ||:| || | || |||  || ||  |:: |||    | ||| || |: :||||  | |||::|:| | |: ||: | | :| |||| :||::||:::|: ||  | |||    | :||:| :|  || |||||  :|| | :  ||  || :| :  | ||:| ::| ||||| ::| :| || |||:|  : | ||: | || || |  || || ||| | | || : |  ||:|| :: |:|   |||| ||  |: :|| ||:|  ||||   :| | ||| :|| |:|  :|::|||||   : || |||: |: |||    | :| ||: |:||:: |||:||||:||| : :|| ||:  |::| :|| | ||:|| | ||::|::::|:::| | | :|| ||||||      |||  :  |:||:|  ||| || |::||  |  || ::| | ||| |:| |:::||::| | :|:|||  || | ||:|  |:||:| |:| || ||| || :| || | | ||:: |  |:||  |||: || :||:|| :| ::| | ||| |:| | |: |||:  :    | :|:  : ||  |:||||  |||::| ||||:: : | |::|: ||::  |:| : |||:|| |::|: |:| || |:|  | |||| | ||  ||| |  ||   |:| |   || ||||:|| | : :| ::||||: :|:||| | ||::| |  | ||| |:||:  |:|||  |||:|:||:| | || : | ||| ||:| |  |:::|  ||||: :|||: |  ||| ||| :||: : :| |:| | | ::  ||  ||||||||||| |:|| :| | || ::| ||| | | || |:| | ||| | | :|||||| ::||  || || :|| :: |: || :| |||| |::|| |:| |: |:|| ||: ::| |  |:|: ||| ||  :| |:|| | :  |:||   | : :| |::|  ||::|:|:| |||| :| ||  : :: ||:||::  || :::||:| | || | :|||| :||:|   ||::  |||  |||:| :||| | :||: |  |:|  | :  ||:| | : | :|| | :| :| :||| :|:| || || | |:|  || |||:|| ||::| ||||| | ||||| : ||||| ||:| :  | :|   ||:|:: ||||: :  :|:| || ||::  |:  ||: |||||| |:|:  :|  || | || |  :||:| || || |:|::||  ::|  :|||||  |:|  |:||  ||| | :|| ||:  ||| :| | ||| |  ||| |:|:  |:||::||  | | |:| :|| :|| |:|: ||:   | | : |||:|| ||: ||: ||  ::||  ||:  |:: || |: |||| :::|   :||:|||||||:: | :|: ||:|| :   :|||| |: |:||:|||  :|   || ||||: :|| | ||  :|||:: : |||: |  ||:  :||:: |:|:| | |  || | | | | ||| ||   | || :|||: |:|| |:| || || || |  | |:| :|| ||: : ||| | || |  |||:  ::|:    |||: |: ||| :|| |  || | |:| | | : ||||||||  |: ||:|  ||:| || |  ||| ||| |:|| | :|:|:|:||: | |   || | |:||| |||   ||| || ::|| :| | | :|  :| | |||| |  ||  :| ||:| :| | |:|| :| |:  |:||| |||:|||: | | : | |  | | ||:||| | |:| |:: |||   | | ||::::||: |||:|: | ||: | |  :||:|| : || :| ||:|||    | ||||| ||||   ||::| : : |||| :||:||  |: |||:| |   :|| :|:| :|: :::|| :|:|:: ||  ||| | ::| :||  |:|| ||| |||   ||:| | ||::||:||    ||| |:|:| | | | ||:|:|:| |::|||  ||: |:: |:| ||||:|: |:|:| |||||||||: |   | ||| | | | || :|  | |||   |:||:||: :||::|:||  :||:|    |||| |::   | ||: || ||||||| | |||| :|: ::||    :|  || | | :: |::| : |||| | :|||: :||   |:| || ||::||  || :| ::|  || |  ||| :|:: | ||  | || |:|  | |:|::::| |:||| || : ||| :|| :| | |::|: |:: | ||| ||| | ::| :|| || :|||| ||: :| :||| |: |:|| | ||  |: |:    |: || |:|| ||||:| | |||: |:  ||| |:|||   :||: : :: || ||::| | | | || | | : |||  |: || | | ||  :||  || :|||| |||: || ||||  :|:||:||:| :|  | | ||||| || : |:  :| ||: |::|| ||  |::|:|   ||   | :: ||| | |||| | : || | :| || ||| |::  ||:| :: ||||  :||||:|: : ||:|: |||  | || :: | | :|:| ::| : |:||| || |||:| ||   :: |:|||:||  ||  :|| : ||  ||||: ||:|   |||  :|| :|||:|: | :||::| ||||| |||| | || :||| ||   ::|:||| | |   | ||||  ||  || ||  ||   |: ::| : || |||||||::|    || ||:|| | : :||:|  :||  :: ||  | :|  |: :| || |:||  || | : || :||||||  ||  || || |:|:|| |: :| ||| |||:|  | : || ||| |: ||:|:|| |||  |||: ||  | |||||| ||:| : :: || |:  |:| :|| |||||:|| | ||| |: ||  :|: :|:||| || :||:  |:  :|| ||  || :|||||:|: | :     ||: || ||:| | ||| :|: ||: | | | ||| : |  ||||:||: ||| ||||: :  :| | | :::|:|| |: | :|||| | | : | | |||||::|| ||| || || ||: ||| ||| || | ||||| || || : |   || |||||:::| |:|| |  |:| |  |  |: ||:| :|| ||  |||| :||||  ||||: | | | | | || || |  |::: |||: |:| :||||   :|||   ||::|:| | | || || || :|  |:| |||| | | |||||  ||||  :|  || | ||||| | | ||| || : || | |||| |:| ||||| | |: ||:|| ||||  |||    | : |:|| | : ||    |||:||| :||  |:  :| |: ||| ||| | | | ||| | ::| | |    | || :||| |:| | | ||  :||::| ||: :|:||   |||:||: | | | ||  :|| ||:||  ||   |  | :|::||: |||| |::|:|| |||  : |||:|:|| |||::|: :|:|| || |:|| :|| :| |  |||  | |:|:| |||:|||| | |:| |||||: ||| |::: || | |  |  || :||||:|:   ||:||   |: || |||   :| ||:  || :|: : | ||| |  || |: || |||::|  | | | ||   |||| ||  |   |||| |:: ||   |:::| :: |||||||  :||: |||:|| :::|  |:| :: ||||: || :|| |::  ||||| |||| |||| ::|: |:|||  :|||   :|  || :|| ||   | | |: : | || | || : :|||: || :|: | |::|| :| ::||  :|| | ||: ||  |  |||| ||| | |||  :||:: | ||||||  :|:| |:|:||  |||   : || |||  | |: |  ||: | : ||||||| :|| | ::| ||  ||:||  | :||::| ||| |:||  |:|| :|  ::||| | |   ||::||| ::|| | :: | | | |  : ||| || |:| | | ||| |||  || || |:||||| :  ||  |:   | || || |:||||  | |||  :: ||| :||   |||     :|| :||||::  :|:|:| | |:|  |||:||: ||: ::||  | || || |||: :| | ||:|||||| |  ::|   ||:| | ||| ||| : | |::|:  |||  ||||| ||:: :: | || ::|  || :| || |:| |||| :|| |||:|||||::| | || | :|||||: || |:||:| || ||| ||||| ::|: ||| |||  :: ||:|| ::: : :|| || | |||::|:|: | |: | :|| ::|:::|||| :| ||:|:||  | ||    | |:| :| ||:   ||| || ||:  :  |::|::|| :|::|:  :||:| : |::|: |:|   |||:| |||:| ||   | :| ||  || |:|    |::|: | | :| : |:|: |||| ::| | ||||:   |:|| |:|||| |   || |:| | | | | |::|:| :|| |: |  || |||||| ||| || ||  |:| :| | |||  | ||| | : |||| ||||||: : :|  :||||| ::| ||| : |: |||  | | ||| | | :|: |||||  :|||  |  |||:  |||    : ||| ||||  |:|| ||| ::|:| | | ||:||| |||: ||: ::||: |||||:| ::: : |::| ||  |:|||: | ||:::: |||  |  |: | ||:  :| ||| |||:  | |::| |||:| |  |: ||: | ::| |||   |:  |:     || : |:|::|:| |||| :  | : :| ||||  |||:| || | | ::| |  |:| | :|:| | |||| :| |: |:||:|| ||||  ||| : | ||  ||||:| |:|||  |:||::||  :| | ||: ||:|  || ||| | |:|  |||  || ||: ||::||||| | | | |  :|: :|||| |:| |  ||:|:||: || || | :|: : | ||  || || | :|:: || |: | |   |||| | :||| |:| ||  ||| ||  || |::::|| || |||| :   |||  :|::|  |::||  |||| |::: |||:| |:|::: | |  ||:  | ||:  :  |: || :|||: ||:|  |||||||| |||: || |:|  | | |  ||: | :| ::::  ||||||  :|:|:  :|: |:| ||  |  ::|:||: ||| | |: || || || | | ||:||  || |||| || ||| || | || | ||  | ||:|:|: |  || | || |:| :|:||| ||  ||   |:: :| |||    || : ||:|||||:|  | ||| |:| ||  ||| : ||| ::|:|:: ||:||| ::||:|||| :  | :||: |||:||  |||| || :| || :| :::|||:| |:|::||  |: |  |||| ||| | ||:|::|   || |||  | ||||  :|: :||:| | ||:| |:|:|:|:|| :||:  |:| || ::|  || |:| :|| ||:|::| | |:|| |:| :|: | || |::|::: | ::|||| | |::|||  :|  | | ||| |||  |  |||:| ||: |  ||: :||  |:  :|||||||:|  | :| : |:|||| ||| ||:| | | ::| :| :  | |: :|:||:||  :|| |  |:| ||:||:|| |:| | | ||| |||::  | :|:|:| || |  | ||||| :||:|:  | :|||  | ||  |  :||| |  |:::| | |:|| |:|:||:| :|| |  | | ::|||  |||  | |   |||| | |: | :| |||  |:::|:||| || |:  :| || :||  | :||  ||  || |||: | || |:|| ||:||  | |   :|| | | | ||:|| |||:|  |:||    ||| | | :| ::||||:|: ||  : |:| ||:| |::|: |||  |||| || |:|| |: ||||||  | :| :  |:|||| |  ||:|::| |  :|:|  | |||:||: |::|| | | ||   |   ||| |||| |:|:||  :|: :::| : :::|| |||: | :| :|| |||::| |: || |:|:| ||:: ||| | ::||| || | :|:||  |:|   :|:|:| |||  :||  | ||| :| |:| | |  | ||||| ||| :|||:  ||::: |  | |:|| ::|: :| | | | |  |:|||||| | :| ||| ||| ||    || | |||:|:|||::||  |  |: | | :|:||| :| || |  :|::||: :||:||:   | |||:| |:  :|:| : : ::|||:|| : || :|| :::|:|| |:| || ||:: :| | | |: |:| || ||||:  |:| | :|| | :|| |: :|:| :||| |||: ||| ::|:|  ||| ||| ||  ||  | ||| |||| ||    || |:|:: :||   | ||| :||| ::||| |:  |  : | ||: || ||:| | :| :| : |:  : ||:||  || : ::| ||| |  :  |  |||| |:|  ||| || | ||  || :| |: |::|| ||| | ||  || ::::|| | :|   | |||   :|::| |  ||||| | || ||| || ::|| : |   |  | | :| ||||| ||:| | ||| :|:|: || : | |||||:| :|:| :|| ||| |  |:  |:|:|| ||::| ||:||: ||  |||:  ||:   ||::|  :|:|:|||: ||  |: |:|:||   || |  |||||:| |||  :||  | |  ||  :||||  :   | || | | |||||||  | | ||  ||:| :   ||| || || | ||||::| || | || |:||: :|::: |::: :||  |||| |:||:| |||: ::| :|| || :|: ||  |:|   ||| |  |||  ||||:|| :|:|: :|| |:|: |: |  : |:| || ||   :||: ||:| |: || | || |:|||  | |   ::||||||| |||: | ||  || | : | ||:::::|| : |: ||||  |||| ||| :|| ||:|| ||  ||| | :|::| :||  : ||:|::|| :| : || :|:| |     |: |||:||:: :|| ||:  ||||  ||  :| ::| || | ||:| :|||| | :| | ||: :|||||   ||:| | :||:||| :||| |||    | :  |:| | | || |:| || :|||  :|| ||| |  ||| || |:::| :| || : ||:| |||||| ||| : |||| | |: || :||:| |:|| |::|:  |:|||| | |||||| :|  |:  | ||| |::|  |  ||:|||||  :| | | :| || : | ||| :| | |||  |:| :  ||  | :|  ||||  :|||| || |  ||: |||| || ||: |:|:  ||:  |:| | ||| |  |||| ::||   | |||  |||| :|::|| | |:||| ::|:| || :   :|| | ||:|| | :|| | | | ||:||  ||:||:| |  | ||  || |  |::|::|:|  | |||  |  || |||  :|:: |||:|| || | :||   :| ::||   | || ||||  |: | | :|||  |  |||| | |||:||::|:||::|  ::||  | |::||   | ||:| ||||: |:::| ||||: ||| ||  :::| || |||  ||   |||| ||:|| :  |||  ::||: :| ::||||||  |: ||| |: : :|| | |:||| ||:|   | |||: ::||  : :| ||:|  || || || : |||:|| | ||:|| || | |||| |||| |: | ||:|  |  :|  |:| | |  ||| || |  |||  ||     | | | ::||||| |: |:|| |||: :|| :||:||||| : :   ||:|:|| || || :||:| |: |||| :|::::||  | ||| :| | |: | |||  | :  :|:||||  |||: ||  :| || ||| ||   | ||  |:||| || |:| ||| |||   |:|| | :|:: :| |||: : : : || ::| || ||    |:| | |:| :||||  |  | | | ||:| || |||| | :| || ||  |  |:||:::  | |:|:| || :||  | |:|||:  | ||: || |::|| :| || | |:|:||||||::| ||   ||: |  | |:: | ||||:   ||:|  |||:|  || |  || || ||:| | ||:: | |: |:::::| || | |   |::|||| |::::|| | |:|  || | ||| ||:|:| :||:|| | : |||  ||:|| ::| ||| |||  :|||:  :|  |:|  :|: ||    :|| ||| || :: |:|| | | | :|:||:   | :||| ||: |:|| ::| |: || | | ||  ||||:| |  |  |:|  ||| | ||:| ::|| || |::| | || ||:|: :  ||| |||:  || |:|||  || |:|  ||||:  |:|:||: :| || |:| |||:|  :|||| | |:|| || :|  |||:|| ||| |   :|: | |||: |:|:| :|:| : | || |:  :||| :|| |: |:|||| | |::| :| | |||| | : | :|| |:| | | |||| | |:||  ::| |::|  ||:|||| |||  ||:| ||     || :|| |  |  |:| | :   |:||||  :|   |:| |:||||||   | ||| ::||: ||| || | ||| |: | ||:|:| |||| ::| :| || |:|  ||:||  |:|||| ::|: |: |:| ||| |::| ||||:| || |  |||   :||::|| :  |||| : |  |: || ||||||||: | |:|| ::|:| | ||:|:  |||  | |||| :  ||| |:|:  |:|  |:|:|:|| ::||:   |||:| ||  ::| :|| || ||   || | :|||||: ||:| |: :| :|:|    ||||::||  :|| |||:|| |  : |||:|||  |  ||:| |: | || ||| | :| |: || | |  :||| :|::||| | ||: ::  | :| |: |||||| |:|:||| | || ||||:|| :|   : | ::|| | :: |: || ::||||  || |||: |: | |:||| :||  |:||| ::| ||:|| :::| |:|   | | ||| :|| | :| || |||| || || ||:| |::|||| :||: :|  || :| | :|||| :||    || |:|| ||:  |||: | ||  ::| ||:| || | | |: || || |  |:||||  ::|||||||  ||| ::| |||  || || :|| |  ||  ||||| ::|||  :|: || | |  :| | ||| :||  |   | ||  : |  ||: |||| ||::|:|:|  :|| |  |  |||  ||:|    ||||:||||  ||||:  | ||| ||:|  ||:|| :||| :|||| :|:: |: :||| |  |  | ||| |  |:: :|:| |  ||| :|| ||:||||  : | :|||| :|| | |  ||| || |||:::|    |:|||   :||| || ||| | | | |||||| ||  | ||||| |:|| :|| |:|: ||:||| |  || |:| :|||::||   |: :  |||| | | ||: |||:::| : | |:|:|:   |:|| ||  ||: | |  |||| |:||  :||     |:||:|: :|:|| |::|:  ||:|    | ||||||   | :|||:|| : |:  |||   || | | ||:| ||||:| |:||  |:|: || | |||:| ||:|| ||  | |  |::| |  |||:| |     |:| || :|| |:| ||: ::|   ::| |:| |::|: : ||:| ||| |  || ||: ||  ||::||| | :   |::  |||| : |||||||::||  :|  :|| :||  | :|  |||  |:||| || :|:||| |    :|:|:::|: | || :  |::|||| |: ||| : :   :  :|  |:::||| |:|| |:| :||| | |||  : |||  | : || ||::|: ||:||  |:| || : ||| | ||| || | | ::|:|: ::|||||:: ||::|| |  ||  | ||:|  ||:  || : : | |||: ||:| |:|:| ::: |  :: ||| ||:|:::| || |  : |  |:|::| |||  |||: || | |  |:||  || |: | |:|:||:||:  || | |  |:|:| :|| || |||  |  | :|  |:|| :||| | |: |||||| : ||:||| || ||  ::||: || |||  :||| | | |  |:|   | | || :||:: |:: | :|  ::| || |||| :| ||| :|:  | ||| :| |: : ||| |  |  |   || | :|  | |||||:|||:|:: ::| :  :|| || | | |:|  |:|||  | |  ||:|| ||: :||  || :| | |||| : | ||  ||| | '
    '-GG-CA-CATAGTGCAGTTACCGTCGACAC--CC-C-C--A-TCCGAAAGACTGACCACAGATA-GTGCACTCCCCCTAAACT-CGC-CGTTCAACAACATCTATATTATC--CAAAAATTTAAAGTTAAAACAT-ACCATAAT-TCTA-TTCACAGA--ATAAATCCTTTCACGAGACGTCC-AC-ACGCGTTT-T-CTCATCCCAGCCCTACATACCCGGTAAAATACCAACCTCTCAAATTCATAAGCGCACC-CTCATCCCTCCAAGGGCTGCCTG-CACGGC-AAGCGTCCACTCCATTCCCTGGTCTCAAGAAGTCAGTCACTAATAACTCAATG-ATCGCAATCGTTAAACGCTGTCCGTAAAC-CCGCCTC-CT---A-A-AAT-CTTCCCATACGT-CTCATAGCGCA-C--CTCAACCTGC-AATATCAGC-ACAAAGGAGCCGACGTC-AAT-CAAAAAATCAACATGTTTCA-TAAC-TCT-C-CCTT-CAAACTCAT-TCCTGTA-CACCTTCTATCA-CA-A-GCGCCC-TC-CACTTA-T-CTACTCA-TAATCACTAACCA-TCGTAAT-T-CCACTGG-GC-ATGGCG-CATATA-TTCGGAT-CGCCGATA-CA-CT-CC-ACAAGTCCCATTAAAC-GCAAAGACTT-AG-GCTCAATTATAATAAGACCGTCTCCCCAAGGGAGATTCACAAATTTGCTTTA-CCT-G-AT-AGTTTC-CCTCACCTAAATC--ACCCCC-GAG-CTGAA-GCCAT-CTGAGTTTCAAAGTTAATTGCTAGTATATAAAGCTACTACAGTTCGA-AAAAGT-CC-TTTANTCAA-CTTACAA-AACTA-AACCATAATATTAATGTAC-TG-T--TACCCAACATACATATCAGCGTTCCC--CGC-A-GAC-CTATCTTAATCTACATCTCTCATAACCACAGGTCATACTTC-AATGAGATCACAGTAT-AGTAAAACGTTTACTTCATCTCTAGTAAGATATATAATGCATCCACGAATGCCCATTGTACTAC-CAGA-CAACATAAATTCCGGAAGGGGATGATTG-TCCCTTTCTACGCTTTCCTCCCG-CTCTCAACCATCACAAATTAACTTAATTACC-CCACAAGCGAAACTCGAGC-CCGGAATACCGAAAGG-CATCAAGAC-CAAGCGC-GTAA-GGAC-CTCGGCCGA-AT-C-TCGTTCAGTCA-A-CATACCCTAATACA-AA-C-TCCCAAATACATCTGTCATGCCTCTAGTCACCGTCTGAATCCAATAATTTTCCAGACTGGCTCAACCCCCCCCTCCCCGTCAATAGTT-C-ATAGGCCACTC-CCCAT-CTGTCGTAATAAACCAT-CAAAGGAACAACATTC-TCTTAGACCCTATATGAAATCACATAACCCCTCCGTACGCCGCCATAACACTATAACCCGACTTTGGTAGACT-CCGCAAT-CAGT--CTGTATATAATCAAATCCACACTCCCTCTGA-ACC-AAACAGAACTATATTCACTA--ACATAATATAAATAAATATTAACAAGCCCTAC-AGTGTAGTCCACTACTTTCTCAGAAGGGTCTAACAGCTTACAGGACTAGACC-CAATTATTCATAAACTCTC-GC-CCC-G-CAATTTCTTCTCGAGATAG-TATAACGAGAAGAATAGAAACCACTATAAAAAGCTAATGAAACAACAAGTAACGTCCAAANACCTGTGACCCCACAACGCTGTTTTTTTCATCCCACCAAACTCTAAGTCGTCAGAAGTCTACCGTCAAAACAACCAGTCTAGC-CCAACATTACTTATAAACCCGACAATTTCACACTGACACTCAGAAG--CTTATACCGCCCCCACG-TGTTCT-AGC--CTAAAGCTACCCCAC-CCACCCCTTAACCCATCAACGTA-CGACCAAACCCATAT-CCTAAACATT-AGGAAC-GC-TACCCGCAGAATGTCCCTATTTA-GGAGCAACATTCAAGACAAATTGCCACCCGACTCACTCGTAGAAAACCAAAAGCACCATACCAAATC-AGAAACTGTGTTTCCATGCATCTCAAACCCTATTTACAAGGGACTTCATAAGAC-CAAAGATCTTTCATCTCGACCCCAAAGCTTGCAAAATC--TC-G-CC-ACCGCTCCCATATCCCAAATGGATTCCTTGTACAGACATACGACCGTTATC-AA-CAATTTTTATGCTT-TCT-TAG-CCT-CATTCACACTTGATAACTGCATTCTCATACCTACACGCACTCCTACGGCGTT-ACCAAGCTCG-AG-CGGAACTCGATACCAA-ATCCTCATACCAGTGCTACTATGGC-GATAAGCCCACACCGACCGTCTTTGAAACTAGGTACAACAATTTACTGTTCTCTGATAGCCCTTCCATAAGATTTACCTACACAACTGGTAAAAACGCACCACGGTCCATACCACGCGGTCCGTGCAAGGAAACCC-TA-CATACAG-AC-TCCGAATAAATT-AAAC--TCTCGTACATTT-ACACGGTTCGACTGTAATGTTATACTACGCTATTCAACACTCATTAGAAACACTACCCAGCGAAGCT-GAGAGACCTCTTCTCATATGCTACCCGGAATCTTTACGTACCATACATATTTTATCAGCTCCTTGCAAATGCAACAC-AACGAAGAATC-C-GTC--CCAGTTCATGACTACTAAATAGTCCTC-ACTTTGCATAAGTTCAACTTATGATCTC-CAAGTATTCTTCACCCGATTGAATATCAACCATCTCAAGTATAAACCGTTTAACCCTTAACCCCCACCACGGCCGACACAGATACTACAGAA-GGCACATAAACACCTTATCGTTCATTTCACTGCCACAGAAACAAATATCGAATACAATTATGGACC--ACCGTAAAACAAAACAACATAAGCCAAATGTTCCACCATCATGAAAAAATACTAAC-AGCCTCA-CTACACAGGCCAGCAGGACTCCCCAATTCCTCATATAACTTAAAACT-GCCATGCGTGACAAAATACCTCATCCACAA-ACCCACCAAATATAGATTCCGTATTTCCATTCAGTCAGAGTTA-CTCAAACCCAGCAAACACCAGCCCTAATTAAAAAAAATAGCCGACGCGCACATCCAACCCGGATGCTCCA-GTGATCACGTAGTTTCTAC--ACAAGTGGC-GAACTTACGTAACCTACCCATACTAAAAAATAAGAAGC-AC-AACCTCCTGTA-TCCG-CAACAGTAT-TGTTTGCACTGGGCACTNAATACTTCGA-AGATAACATACAAT-ACTCGGACATTCTAGC--TCCCCATCCGCGATTACTAATCGCTCTTCTATATTCTAACTGG-AACCCGTCTGTACCCGCCTATTCACGGACT-ACCCCATACGTCGGGTAT---AGCATCA-GC-TAT-TC-CTTCCAATAATTCCACCGCGCTCT-CCTTAGGCCAC-CGCTGGTCAACACTCT-TG-CCCCACAACTTACAAGGACCATGGCCTAAAAGC-CCGAAATG-CGCGACTCCTAGCGTAGGAGGTACGAGTCGAAACTCAATCATCT-AGCGCCCAGCTTTCT-CATGAGCTCAAATGTTACGATTCACCTCACATTACA-GATA-GAAC-ACACTGCCACCTTAGAAGCTTCC-AACGCTCAGACAGCTATTCAGCTTTCAACTA-TGCT-GCTGCACCATGTCTACC-GGAAACCTTA-T-A-CT-C-AAGCACGA-TTCCGTAATCCG---C-G-TACAATCACCCCGAC-CATCCTC-C-CAGTTACTGCCCTGTCTACGGAACTGTTCTTTGCTAATAG-C-TCCCAATGGGTATTATGTCAA-TCTTCAGTAAACTTACTC-A-CACT-TGT---TTAC-C--TCCCCAAATGCAACA-ATCC-TTCACA--ACA-CGC--ATCTCTCGA-CCAATCGTAGACCAAGTACCTA-CTCCC---AGTG-TTA-CCCTTGAACCTTCCATTGTGGTCCT-GACCA-AACACG-CCAC-CCCA-CAAC-CCAACGCATGCCATAATTGAGACGCC--CAATGAC-AATATCATACCATC-CACTCAAATTA-CAGAAGTAG-TCACA-A-TCTAG-CAA--GCTGGGTTGGCAAGTGCCACCCTG-CGCAAAAACACATGTAAGTTTGTATCCTATCCCACATTTCATCTGACCTCACACCAAATGTACTCAACATCTTTTAAGACATTCTGCCCCCCGTACAACCCTATCAAGCGGTATATACAGAATGCAAGCC-TTCAATCTTATCTTCCGAACGCACAACAACCTTCCGTATCTACC-ACCAGTCACCAGTGTATAGAACCCGATTATAGGTAC-GCCTCACTACATTATGCTCCCCCATCAGCCGCCCCCAACGATTCTCACCCCACGTTTTAAACGTCTT-AACAAAGCTTCCATCACCACCGTACTTATCTCCACCTAATACGGCCCCTCAGACACGACGCCAGACACTCTCACGGGCAATTGTCTCAGCATTATGC-TCTCTTATACCGTCTCTCACATCGAATTAATCTAGTTCACC-ACTCTCTTCGG-ATCTTGAAAATCGACTAAG-TAAGAACACA-GACT-TA--ACCGA-AATCGCCTTGGTGCAGCAGCTTTCCTNCGCAGAAATTGCTC-TTCCGTAGCACTAACGAAGTCCTTTTCTCCATCTCCATG-CATCTCTGACAAACAAAACACCATGACCC-CCCAAATACCCCCCTAAAACTGATTCCCACCCCA-C-ACAACATACTGAATCGGGAAACTCACACCACGCCTAATATCCGCCTAAGCATTCTGTTGTA-CTCAGTAGGACTGTGTTAGCAACAATAGAAGCGACAATATACGACCCCGCTCGAATTGC-CA-TA-CCTTAC-CTAACCACCTAGATCACCAATCCA-GCC--TACCACCAAGCTCGTAC--TCCCGAATACTAGTTTCCTATAGGGTGAAATATGTC-CCACCATCAATACAGCTCGAATTTGTGCCTGTCGTCACACTTGCTATGCTTTA-TGCTTAACCATCAGTTCG-CAGACTTTCACGACATAACTTCACTTACAGAGAAAAACGTCCCTGAACAGT-C-TCGATAAGCACAAGGACACCCGCCGCCCT-C-AAAACCCCCGTC-GCATTATAGTCCCCTGCTTCTTATACACATGGCCCCCCTCCACCTTCATAGTCTAACNCCGCCCAAACATCTTCTGTCTTACCCTGGTACGACAATCCTAGTTACTGAAAACTCCACCC-CTACATAAATTGCATC-CCCGCGCCCACAAACACGCTCCTATACAAAGCA-AAATACACGATGGGCCTATTAGC-TCTCTACTCTCTTCAACAAAGACCGATTACCCCAATCAAGC-CTCCACTCCTGGAGGCCCAAAG---ACT-CTTGAAATCTGATAGTATAAATTACCA-ATTCTACAACCTTTTTCAATTTGACCCTCA-TCCATCTAAAGCTATCCACGTTCAGTGACACCCGCCCAGCTTCCTGATAC-CCTTACCGCATCTTGTTCACGTAGGTTCCTTATAGGTTCCGGACTGCAA-CGCTTTTTTCTCATCCATGCACGCCACC-CATCGGCCACTGACAATTGTGATACGCTCCATATAC-ATCAATACTAAACGACACCTCCGTAAACCTA-CCC--CTAACAAAAT-TC--CATC-CACA--CAATAATCGTC-TCGCTATTTC-A-A-CCCCTTCACCC-CGCACGTTAGACCCTCCCCTGACCTCCTATATTAAGGAAGCACTACTACA-ATGCCTAGACTCCCCCCTGTCACCTTCT-ACTCGCCAATCACCACT--CGTAGCATACAGAA-TCCTCGTTTAACTATTGCTACACCAAGTCCGACA-TCTACACCGATGGCTTTCTCCCCCATTTTAGCAATCCGTACCGCCGATAACA-AACTACGTCACCAGTAGATCAGAC-ACGGGAAACCATCCCATCTATAAGTAACCATCCCACCCCGATTAA-ATAATGCGTTCATACAATGCATCCCAGCGACGACCTGA-ACGTCTACGTTACTACGCTTAAATGCACTTCAGTGCCCGGCA-ACCACTAACACGTGCATCG-CACTATTGTTTTCTGACGCAATTACAGAACGCACTTAGACGTTTTATATCACTACTACC-TATTGCTTACTGCATTAACCTCTTCTTCCGTTCGGTAAATCACACGTTCCATCAA-GCAAAGAAATCTAACTAACACT-C-CATGCC-GGTACCCTTCTTTGCCATGTCCTTACTCCCCAATG-CACCATCTCACCCTTCCCTCCTTAAAAACCGTTGAA-TTTAAGCTTCCCCCG-AACTCCCTAGTCATCCCCGCCA-GCAGCACTGAAAATTACTATTAAACACATAACCGCG--AAAC-GAT-A-CAATC-CGCCAATGTAGCAAACTAGGAACACCTACCCGTTGATGTCTAACCAACCAGCTAGTTGCGATTTCTTCCCCCGAAGAACCTAGAAAGCGAACATTT-ATTTTACAGAAAGTACTTTC-AGCAGCCCAAGAGTACCCTTCTCGTACGAAATCTACAAATAAAATTTCT-TGTACCTAAGAT-TTA-CCCTT-CTGTAGAATAACACCACCCCATAAACCTTTTCCAGAATCGGTATATTTTTAATATTTGCCTATGCCCAAGTCTCCACAAGAAT-CTTCACGAAGCCTCTACCTATAAGTAGTTATATTTACCGACCTATTCCCCCAATACCATCCTATACTACTTCAATTCAAAGCAGTACTTGACCAGTACCC-GACAGTCTTGCCGTGAATCATCGATG-CTCT-ATACCCTCACAGTCAGAACCATACTG-GCCATCAGACCC-GAT-AAATCGCATCATCCGCGTAACCCACCAATGGCAATGTACATA-TC-T-CCGTTCGATCCT-C-ACCTAA-C-GAA-A-TTATTTCCATTGAATAGA-GAGGATATATGTCCATACAACCCCTATGAAAGATTACAGACGCATGTTGTACCAACTTACCTGCTTGCACAACCACACTCAAACCCAGTCTCATGTTCCACAAGTTTACCTCAGACCTTGCTATAGGATTCATTGCC-AAACAACGCACATACTATAT-TCTTTT-TCAAACTCTTATCCATTCACTTCATATTACTCCAATGNC-TCCTTTG-ACAGACAGGTCCACTTACATTCTTATCAACCCT-GTTGT-AGCTGACGAACTAGGCTTACTC-CGAAATCTACCCAATCATATAACCCTTC--CCCTGGTGTTACATACTTAAC-CTACAGAAACAAATCTTCAACGCGTTGGGCTACCCGCCAATTTCTGCGAACGTACTGCATAAAACCTTACCGCGGAA--CAAT-TACCGACACAATTTAACAGCCGGCCTC-CGTTAAT-ACAACCTTATCAATAGTAA--GCGTTACCAAGATTCTCTTT-CCT---CTA-A-CCCAAACACA-CTT-TGC-TTAGT-AATCTTTTCGCCCACCAATTATAC-CTTATTCCCAAACAA-TGGGA-ACAATAAGC-ACCTAGCTCAGTTAACT-AAAC--CAGCTTAACTAA-CCTCT-TCAGTACTGCTAAACC-CTACGACAAACCATAC-T-TGCACAACCCACTGACC-GTTCTCATC-GATCAGTTTC-TGT-CCTCCGAATACCANCAT-TTCTGGACATTACCTACATGTCATAGGGAACG-GGAACAAGCGCACTTACCCC-CAT-C-TAATC-CCG-CCCATTACAACAGAAGAATTAACTCCCACGCTTGAAGTCATCCGTGCACTCATAGTC-C--CAC-G-CAACC-CATG-GAAACCCCGAACCGCTA-CAACTTGC-CCTGTTCGTACAGTTCCATTACTA-CTATCC-TCTGACTATAAT-ACAATCGTTTCTCACCTTTGATC-CTTCTGTATCTTGCAGACAC--CCCCCTCAAAACGA-AACAAATCTGGAGCATATACCTGAAACACCGGACATCTTAAAGCCTCGCA-AAAC---CAAAACACTCAATT-TATTCATCCA-CATTC-CTGATAG-AAACATTAGAATGAGACTCAACATCGCCAGC-TTC-GTGGCTACG-AACCCACCCCGCACTCAAAAACTAACTTATTGAATTCAATCTACATCTACTTTTGAGGCTACTCGATCGGCCAGATAAATTATTACGC-ACCATTCGCACTAAAC-CAGCCTACACCCCAATTCGCTATTTAAAAAGTTCAATTGCACC--TCCGTGATTTTCACCTTCCCCCTATCAGTCAATTCC-ACATGAGTCTCAAATCCCCGATCTACACGCAACCTATTACCCCCTACTACTGCATTTGTTAAGCCGTTGTATCCAACAGTAAATTTACCGTAACTACTTCCTCTGCCCAAGCAG-CCCTACAAAAACTTATACACCACCCCTTTAGTCGTGTAACTCGGTCTCCACCCGTCA-AAACTACATACAAATGACGAGTAGGTAG-CTATATCATAGTATAACATTCCTTACATATACGTTTCAACAAGTCCACATCCTTATCACAAAGTACTATAAATTCCAGAAACTAAAGGAACTCACGATCCCGTAGTCTGAGTAACAAACCGACTTTTCACGACCG-TA-CTAGCCATTCTATCTGTACATGCTCCATAT-GT-GTCTTTGC-CTG-GAACAGCCAATTGCTCCCCCCGATTCGTCATAGC-GCAACACAACA-ATACTAC-GAT-GTACCCCAACTGCTCACGGTCCTATCA-ACACTATTTAGTGC-ATCTACCGCTTCTACAAAGACTAACATCA-ACAACATACCCAA-TACTTACT-CTCAAAAT-AGTC--CGCTTCCATA-CTACTCAACCAAAACATACTGACTTTA-ACCTCCCCAA-CGCTAAAAAGCCCCCAATAGACCACATTTTCTCTCCCAA-ATCGTCCNAA-CCGAATTGAACACATTTCATCAAACCTCTA-CCGGGATCACAAATACGGTATTAACCGC-CGGTCCCTCAACCCTCTACCCTA-ATAGAACAATAATCGCCCTAAGGAC-CTCATAGTTCAACCACATTGGTCCACACTACTCACATTAATTTCG-CACCTTCTACAAGTCTTGAATC-TGAGCTCGACCT-CAGACTT-TATGCAGTGTCC-T-CGACTG-CGTTTATCACGTCACCTCTGCACCCAGGGGAACC-GCT-ACACCTATTCGGACCCCCAAACCTCAAGC-TCCAACTCTCACACCGTCACTGT-GCG-ATAGCGCTA-AGCACACGACAATACAGTACTGTCTCCTATATGTGATCATTATACTAGCCAACAAAAAATCCACCCTACATTTCGACTGATCA-TATCATTATACCCTATACTCCGTCAACTAAAGTAACACTGCTGTTCGACAACACGAAT-TAT-ACAGTCCGTTTCTA-TCACCAATATCCCAGACCACTGCGCCAAAACACAGGATACTGACCTTGGCAA---AAGGAGTTT-CCCAGCTACG-AACTAACACAAGTA-AGCAGCG--ATGGCCTCTTCACACTCCAACAACGCAACTACTCCATAGCCCC-CAT-TCCGT-ACATCTAAAGATTACCAGGACCCCACTAAACCAT-A-GTCCCATTCACAG-GCA--AGCCACCTG-CTT-CCG-T-T-C-CCC-CAAACTC-G-CCCCTTAGTACC-TGCGCCTATTATTCCTGT-TGAGTCGTC---TCCTTATAAGCTGCTGTCTAGGTTAAAAGGCAC-A--G-AAAACAAATTCA-CGAGACACACTGCTACAATTACCTGATAATACAAACTTAC-CTTTCATGCTC-CCCCTCGGTCTTCAGAACTAAATCTAGTAC-ACCCTGTCGCAAATTCGCTTAGCCAGCCGCTTAGCACG--GCTTAC---GGCTC-CCACATTTATATATATTTTAAAAATCCC-T--CCAAAAATATAAAAG-GA-CGTTTATTGAGAATAGCCCCGTCACC-CTTAGTG--GAAAC-AGTACTTAAC-CAACTCACGTGCTAGTACCTACTTATGTATTGTC-ACCTTAACTAATAACCTAG-CATACAACGTCACTA-CGACC--CTC-CTAGATTCAC--CCACGAGCACT-CG-C-AACGCTCATAGCCCAAAGC-CAAACAAC-AAAGCCTCT-C-CGGTGTAAAGTGCAT-CGAATCACATAAACAATG-ACCCAGCCTATA-CACTACACCGCCCGTCCATAACC-T-CTCG-CCCACT-GCCTTACTACTAAGACAAC-CACACCAGCACCGGGGTACGAGAGCGATACCTGACGTCCAAACAATATAAACCATAATAGGTTCGCTGGCCCTATGCCATACT-GCACTACTACAATGTCTCTCAAGAAATTAATCAGTTCCAAGCAGATTTGCNCCAGTTACTACCACTACCTGCCACGAATAGTAGAGAGGATACAAACATTTCCGCAAGACTGAAATGGCCAAAGCCTTCAGTC--CCCATACGATATCTCTGCTGTACACACGAAATGTGGAACAATCATTTACCTCACAGCGAAATTCTCATACCAGCAAAACAACCAGCCTTTCCAAGAACCTACAGCTACCGGCTCACTAAGAACATTGTGGAGCTTCAAGAG-AC-AGCAGAGTCAATTGTAAGTGAATAAGCTCTACGCAATGTGTAAAAATTCAACGTCTATCCCAATTGT-CGTA-CATACATATTTTTCTATACATTTCAAGAAATCACCCCACGTCCGTCGCGCCCCGGATCGTAACCGCCCATAAATCAAGTTCAAGGTCTACACAAGCTCCACGAGATAATCATGGAGGTAAATTTCGCACGACTCAATTAAAAACATACCTTTTGCCTGT-CGTT-TTTCT-TCAAGTATAACCAATTTAATG--TAA-ATCAC-C-C-A-ATTTGGCACCACACGCCTC-CCCCACGATAGACAATGCCGCNTCATGATATTATCCT-GCACTCATCACATATTTTTTGATACTCAA-AATCCTAAACATCGTCCCACCCTTTATTATTGTCACAAATGTCCACCCCCGTATGCTCCTTCCACTTCATAAACCTACACGGCGGTCTGCTTTAACGCAATCACACAATGAGCCGATAATAACCATTCTCATTCACCC-ATTTCT-G-CTGGATCAACCCATCTTTCGCAATTAC-CCC-CCGATTTGTGTGTGCATTGATTCTCCCAGCCATCCTATCTCAACCATCCCCCCTCTCGTATGCAAAAACCTT--CTTCACGATTG--CACACTAC-AGATTTA-TCCTAACCAAGGCGC-CGCGATCACACAACAATCGCTCA--CCCGAGTTTC-CAAACTC-TACCC-CATCCGACC--GC-ATAAGCTTACTCAC-TTT-A-CGC--AAT-CTC-AAGCACTTCCATGTACTGC-ACTAAACTTCCCTCACACGAATAATCCACTACCG-TAACAGG-TCATTACAAGAAAGTCATAAC-CAGAGCTGTTAAATTATACT-CTGACCCTTAC-CTAAAACTATTACC-TATAAGTCACCCTACCCGGATCAAGTTGGTTTACATAAATACCAACGCATAAACACTGAACAAGATACAAAAGTATCTAAGTATAGTGAATTACAAGCAT-AGATCACAAATTACGACACCCTACCCTCCAAATCCGAACCCCCGCGGACTGAAGAGATACAGTTTTTCATTACCCTGTGCACCCTAC-GC-C-T-GAACAG-CCTGATATCA-ACA-CATTATGTC-CTTCT-CCTATGTCTTTCCCACACCAACTAGACCCCAACCACCGGAAAGGGAATCGGTGAATAATCACTAACACTACACAACAAT-GCTGTTACATTAATATCACTTCTATCCACATTAACAATCGCCACTCAACGCCTTATCAAACGTACAAGACAAC-AAACAAGAACAAGTACAATACATTATCC-CGCTCGACGAAGACGCCCCGTCAGAGTCCATAAAGTTAGACCACCAAACTTAAAACAACGGCAACCCTGCCAAAGC-CCACAGATTACGCAATTACAC-TACA-CACGTGTCACTTACTACAGACTGTTAAAAGGCAAACCC-CACTTAATCACTCAAA-AAGGTGTAT-ATACTCCACATACT--AAG-GCTTAACATCC-A-TCTTTATCTCTACCA-AAAG-CATAAACGAGCTCCAAAACACTAGACCTGAAAACGGAACACATCGGAG-TAC-TAACCACGTCTAATATTTGTAAC-TC-A-CCCGCATTCAAAACATAC-ACATC-CAACCCGATATCT-TAAAAAAGCCCCTATCCGCTCATC-TCATCTTTCGTTTCTTTAAATGTC-T-CCCCCT-G-GATCAATATACTACAAATNATAAGATATCCATCAA-ACAGGGAACACTAGATCCTCCCAGCTATTCACCCACTTCCCTGATATTAGCCGATATCCTCATCTCGCTTGATCTTAGGCGTAGTAACGTTCCCAACAGATAATTCAACACACCATGAT-TTTCCCAAACTTAAAATCAATACTTCCATGGAACACCCTAGAACCTAGTAGATACTGGACTTCCA-TTAGTAAACTAGCACCACTTTCCACCCACTAACTCTAAAGCCGTGATTTTTCTAGAAATACCCCATGGAATTCCGGTTGTACCACCCCAAGACCACAGCAATAC-GTCATCAATCCAAAACGCCACCCCTAATAC-TTCATGCATCGCCACTACATAATGCATTACGATTGATAACATACATCAGGAACCTCTGGACCTAAGACTTCGCCTTACTTAGCAGAGCGCAAATGTC-CCCGTGTTTTTCTCTAAAACTCCCTTAAC-ATCCGCACGTGTTATAGTACTTTCCGACTAATACATCGATTATTCTACCAACGTCGACAA-AGGCGGCCTTTCCGCACGGTTCTTTCATGTTATCTATTGTTCAACGCCCCTTCAATCCTTAATAACATGTCGTAACT-ACTGATTTCAT-ACAA-CAAGTGGCCAAAGAGTAATACCCTCTGCCTAGTGAGATATGAGCCAATACACTTGAACGGGTATATCAGTTTCCAAAGCCACAGCCGTTTACTCGGACTCTTCTCCTCGTGCTCCGCCTTGCTGACCTGCGCTACCAACCCCCTCCTCTTCACACTATTATCATTAT-ACTGCA-AATAACAAGCTAGTGACCGA--CGTTTTACCGCAACTAAC-GTCC-CAATTCTCTCAACTTGATAATCTC-TCTGATATAAATTATGTTTACGTCCGACCGACAACA-CCC-AANCTGGTAATCGAACTTTCTACTA-ACCTGACACTAAATTGACTTA-CATGTCCAGTATCAGCCCTTATACCTTCTCATGCTAATAATCCTTCCAGTCAGATCGCC-C-CT-ACCTCATTGATGACCTTACCTACG-ATCTGGTAGCAAATTTTAAATCCAGCCCCCGGTCACATTATTGCAACACCCCGACCGATGAGAGAATGTA-ACGT-TAACACATCTAGAAACGTTCNCCCCACCACTATCGGTTACAGGAA-AGATAATTATGAGCCCACATACAATAANCCCGGAAACAAGAATGCACTCTCCTCAAAAGACCGCCTG-AACTCC--ACATTAGTTGCATCAGA-CC-AGTCC--A-G-CCTTTGCATACGCGAT-C-TG-CCTGCCACGAACTAATTCAGATAACATAAATCTAT-TAA-AAAATTACC-TAC-TAT-T-AGAACCACTCCGCCCTTTGCTCAC-TCGC-AGAT-C-CCCCCTTTTTTAAATGACATACACCTCCATCTT--TACGTC-CGTACAGTAAAATTCCCCTCCCAAAATGCGACCAGGCCACACCAGAACATAATTCACACATTTATGCGCTCTGACTACACTTCTTCACAACCAAGGACCTGTCTTAGTCCAAATTCATGTCCTCACAA-CTTCTT-ATATAGAAATCA-AACTCCTCCAACTTAACATTCGCCAAACAATTACAAACAGAC-CCGTGC-CTTACGTCTGATAAATATTAC--TCACAGCNCAAACGCTATCA-AACACATTCACGCAAAGCC--G-AC-CCTATG-AAAGTGTTAACCGAACGCCCGAC-CCAGATCAGACATTCCAACGTAAAACTTAAACCTTCCTATCGCCT--TCCCCATAGCC-C-GTATA-GTTAC-TAGAGATTAGTCTCGTCACA-TCATGACTTAATGGCATCCTGCACCCCCACATTCCA-C--CATGAATTACCT-ACAC-GTCTCATTAGTGTAATAACATG-CCTCGTCA-TAACATTTTTGATAGCAAATTTTCA-ATTTGTATG-GCCCACGCATCCATCAAAGATAATTATTATTTCTCATGACTGTATCTTCATAAATCTGTTCACCCTGCCACGGCCTTATC-TCGTATATCCGAAACCA-CTGCA-GCTTGGCCATACTCACTTTATGCAACCCA-TCCTGGGACACGCCCTACGCGTGGCATGCTCCATCAAACGTGCGGCATACCTAC-CTCACCTCCACACCTC-CCTCGCCCGGGACAT-ACTT-GG-GTAAAACCCCCGCGTTTGCCCCCCACA-AC-A-CACC--TG--CCCA-A-A-AACTCAATTACTACT-GAAACTCAGGATGTCGCCCACTTGATAACACAGTATCTTAACACAGTCACCAGAGATAAATCCATAAATATTCCCTCCTCAACCT-TCGAATACCG--GTAAGCTTCTATTCNCCC-GCAGTAGCAAGTTATATACCAACGGCACATTGCGATCCCCCAAAAAACGCCTACTCATATGTGTGT-CC-GC----GTC-CCCAATAAAACC-A-CATG-CTAAAAACC-TTCT-G-TCCGA-ATATAAGTATCTNGGTGTCTTTTATC-GCGAGCTAACACGCCCTGAAGAG-TAACAGAC-ACCA-GAAAACTATCAACACAGGCCCAACCCG-TAGCCTACCTG--CTCAAC-CCGAAAATACA--TTTCA-AATC-T-AGGTCTGCTCCGTAAGCTTCCG-A---CGGGCGC-GTTTGCACG-TTGCTAACCCTACGCCTGAAATTCGCTATCGCACACCGCGTAAAAC-AACACTACGGGCCT-CGCCGCAG-CGCTTACGGAGATCTAACGGTCAGGCTAC-GCCTCTATTACACCCCCACTTCTCACATACGATCACTTCCGC-CTTAGGCACTTACACCCTCTAATTTATCACATGAAATACAGCCATAAGGCCACACTCTTTTCCCTAGCGATCTTCAGCTT-CCAGCGC-CCACCACTCTAACTGCGTCTCGATGTTCCTATTACAGTTA-TATCTTCTACCTTCAATA-TAA-CTGTTCTGAGAACAGAACTCAATGC-TCAGAACCACGAGTC-GAT-CAGTATACT-A-TAATAAATC-TAAATTAAAT-TCCCATGCTGA-CCCA-CTTAGGACGTCCCAATCCATGCAATCAACTTCCTC-AA-CC-CGACAAC-TTTC--CAT-TCTGAGCGCCGCCCGTCCGCTACATAAATCC---CATAGTTTACACATACTCCT-AAC-TA-CTGATCACATCTTTATG-TCACC-AAT--CAGAA--AGACACGTTCGCAT-TAC-TAT-AAAG-ATTCTG-CA-CTGTTC-CCTCCGTGT-TCCAACCTG-TCACG-CATGGAAACACTCTCATACCAAGGCTTAACTAGTATTC--A-TTGTAAAGATA-ATGGCCCAATGCTTTCTCC-TTCAAATTTCGTCTACGTTGCTGATTCC-AC---CC-TCAAACCCGACTTACCATGCTTCTA-T-TACACAACG--AAC-CTTAAATCGGACCAAACGC-CT-CATTCAGT-CCGCGACATCTAAGTTCCATATACTTAACTTGTTGAGACCACAAAATAGAAC-CAGGTCCCCTTTAAATTAACAACCCCTATCCACCGCT-CGATGCCTAATAGTGCCTTAAATG--TTCTTAG-CACCCTATCAGAAC--TAAACGATCTAATTAAAGTATCATC-CCAGAACCAACGAACATTAC-CCCTCACCCATACACCAGAACAT-TCCAGATCCAGAAAAC-ACTCATGAGCGTTAGATGCC-CCACCATAA-TCGCT-AAC-CCACACG-AAATAAATC-CTGT--CACATATCATATCAAAGCCACCTAAGCCCTAGCGTCC-TACCCGTACCCCGGAC-AAAAGCCAAAATCACACTGGCCCTGCTTTCCGAAATTGTCAGTATG-AAAATTTCACCAGTATCACGTTTTGCAGATGCCCGAATATCCCCCAAATTAGTCCAACCGCCTCCCTGAAAC-CAAGCA-ACTACAGTGCATATATCG-GTTTAATCACCCAGACGCTACTACGTAATAGCCAA-GACCTCTCCTTGCTTCGACCCACTAACCCA--CTCTTC-CTTTTAG---CGAATC---GACA-ATTTCCTA-G-CACCAACCCCCTCACAATCT-CACCTGAACAACTACCAACCTTGGACGTAATTAAAAACTTTCGTATCCATACAGGCCAATTGAAAAACAATGCCTCGCACCCGTCACTAGCCAACACTCCATCC-CGCAGGGCCGTAGGCATCACTCACATACTTTTTCCGTGCATTTGCTTTACAGGGCAATCCGCATTCAATTAGCGAGCCCATCAATAT-TAATAG-AATA-TTACGGCTCTGCTATTCCTTCATAG-CTCCCCCAGCTAGCTACTCTC--CCCGACCATTAAA-CAGACCA-GGCTGAAACTTTTGAAACTAAACCAC-AGTAC--GA-GCT-ACAGAAGTGACTCTGAATCGG-CCAGT-GGACTACTTGTTAGTGC--GCCTTCCGGACAACCTAC--TACCCACACATTCCATAGGAGAGGCTTATACTCACCT--ATATATGGAATCCA-ACTCCCACCTCCCTAAACAC-CAAC-GAAATATA-C-TCCC-ACACCATGCACGGAATGATGTTAGTCTGGAATCCCTTCCT-CTACTATCTA-C-AAAACAGTGCAAACTAACAAGTCCCCCCC-CA-ACGC-CTAAATTCAAAAAC-AACAACAGC-CGTTACAAGGGTAATTAACCCAGTTAAAATC-CGCACT--CGAGACAGCCCATTATCATGCTAGACTTCCGCACCAGACAGGACAATCTCGCATACTCCTGTGGGAACTTATCTCTATAGGGGTCCCATTAGCTCTCTCACCCCCAACATTAATAT-CTGAGTTAGCATAAACACTGGTCCTACCCATAAAACCGGACCCCGCNGAGACCAACGTTCCCAAACAATAAACTTCTATAAAGTCCATTTGCATGGAGTTCTTCCAACGCGTTTCCATTTGTAACCAAACTACAAGTTAAAACAGCAATGCAAATA-TGCATGACTCCCT-TC-AATGCACTTGAACGTTGTACCACTACC-A-ATTCCCAACCAC'

>> showalignment(rglobalAlignment)
>> 
>> 
>> hippoPORF_random = hippoPORF(randperm(length(hippoPORF)))

hippoPORF_random =

    'QIGNKMLTFPRIIEAALTTTTLIVQGLTITNILPLTLGQRQIFTNFTI*ALYLHITITSLTISGLLIGLTPVTPVEHILLSTFPLHANNAILPMYTPIISAILLLLMPFLLVSIIQLTIAIQLINRSLNIKLSIINLLTGFSSTAVNIAGISHMALALILTFFPISLVG*LNA*QPMGEIDALIPVTLHPKHTPFVVIQVLSTLLIHALNTIAKGAFQLLFIFRIG'

>> humanPORF_random = humanPORF(randperm(length(humanPORF)))

humanPORF_random =

    'VFRYTLFVREILIISDGTEAQTTDFNIVELPLILYKSITSALMANRIPFANLLL*LSIYLVEIAAIGVPIETLTFYAGR*KPSPITLPLLAVIILISASSTPYAMELLALARAPLG*ITTYLPLIAISVLMLLLGTLTIYTS*YYEFGFPRFTTLHSSLLPSLTYTLSNIPK*QLIIIKLN**TFTNFTSALL*PGTIFPALLPAVLIIHAADNPFPIA*ATGLELGRALLKELNFEIETTYSALLQLYFLVSPILAYPTLTSGANNNGTKFPIYTILGSLLFLALSQKTSTDLNTLLAITLLTLQNTPLIIVLQI'

>> [rscore, rglobalAlignment] = nwalign(hippoPORF_random, humanPORF_random)

rscore =

   -95


rglobalAlignment =

  3×317 char array

    'QIGNKMLTFP-RII-EAA-LTTTTL-IVQ-GLTI-TNILPLTLGQRQIFTNFTI*-ALYL-HITITSLTISGL-L-IG--L-TPVT-P-VE-HILL--STFP--LH--A-NNA-I--LPMYTPIIS-AILLLLMPFLLV--S----I-I-QL-TI-AIQL--INRSL-NI-K--LSII--N--LLTGFSSTAVNIAGISHMAL-A-LIL-TF---FPIS-LVG-*LN-A-*QPMG-EID----A------LI-PVTLHP-------KH-TPFVVIQVL-STL-L-I-H--A--LNT-IA-K-GAFQ---LLFIFRIG'
    ' :   :::    || :::   || : ||:  | :  :|    :::|  |:|: :| ::|| :|:  :: |  | :  |    :|:| | :   ||:  |: |  ::  |   | :  :  | |:|: ::|:||:  | :  |    : : :: |: :  |  :: :| || |  | ||  |   :|:|:|: :  : |    | | ||: :    |||:  :|  |: |  : :: ||:    |      |: |:  :|       :: | | :  :| | | | : :  :  ||| :|    ::|   |:::::| '
    'VFRYTLFVREILIISDGTEAQTTDFNIVELPLILYKSITSALMANRIPFANLLL*LSIYLVEIAAIGVPIETLTFYAGR*KPSPITLPLLAVIILISASSTPYAMELLALARAPLG*ITTYLPLIAISVLMLLLGTLTIYTS*YYEFGFPRFTTLHSSLLPSLTYTLSNIPK*QLIIIKLN**TFTNFTSALL*PGTIFPALLPAVLIIHAADNPFPIA*ATGLELGRALLKELNFEIETTYSALLQLYFLVSPILAYPTLTSGANNNGTKFPIYTILGSLLFLALSQKTSTDLNTLLAITLLTLQNTPLIIVLQI-'

>> showalignment(rglobalAlignment)
>> 
>> 
>> [rscore, rlocalAlignment] = swalign(hippo_random, human_random)

rscore =

   1.7111e+04


rlocalAlignment =

  3×17541 char array

    'CGGTAACCGAAATGCAAACCT-AGACCCAAC-CT-CACAAATTC-CGAGAATCGATACTGTTACCTAT-TCATCC-AGTTTTCTCCGCACTTTCG-ACTGTATAATC-ACAC-CAGCG-GTA-CCCTAC-CG-GACGCTCCTCGACAAATTATCTACTACAAGACTATTAATCCC--ACCGATT-CATCCTTCTAAAAAAC-T--GTGCTCTCTCGCCCAATCTAAAAAC--CTTACCTAACTTTTGAATAAT--ATTCGACAGGTCG-ATAAAAAAAAT-GCTGATTGACT-ATTCTATACAATAGGACATCTCCCAGCCGCGAGTAC-TATCACTGATATC-A-CCAAGTATCTCCTA-C-TACTAAAACAATACTATACATCGGCTCACTAGAACAGAATAATTTAAAAAAATACTCAAG-AT-GT-AC-CGGCTCATCCGACTTGTATCC--CA-ACCTTTTTAC-CATTTACTAATGTCAATGGTGGACG-AA-C-CTCATTGACG-AATGCATCT-AGCG-GACTCATCTCCA-ATAGCGCC-ACAAACTCCATATGGTCCCTTCTCACACATTACACCACAGTAAACGATTCCCATTAGT-AAACAGCCGATAGTGCTAATTCAC-CTTCTACTTCCCCTACATAATTATTATCTCC-CAACCGGTACCTAATTGCTTCTGGCACATGAG-TTCTCGATTACAT-CACATAATTGTCCCTATATGC-GTGTCAAGATCATGCCGGAGATTAGTATGGGGATCGATTAACATCA-ATCAAACCTCAGGAATCAAAACTCCCTA-CA-CTAACTTAG-AGACCTTAGCTTTTTAGTGACA-CAATGCAAAACAACTAAT-ATACCGAATATCCGACCCCACTACCGCTACAATC-TAGACATTAATGCACATACTAGTATCCTACAAATG-AACAGC-AA-ACCCATTG--CTACAACGC-TCACAAGAACT-GTT-TCTAT-AATGTCCGCCTTCCAT-CTTTCGT-GAGGGAAC-ACCATATCTGTAGAAAGATATAAACGCAAAATAGTCGTTTTA-CCCAGAGATACT-CACGACTACGGAGGCCAA-CGC--AA--CCCAC-GCC--TTCAGCGATTGATGCCACTATACGGGTTGC-ATCGGTCACTTA-A-CT-TCACACTTCTAAAGTTCCG-CCGCAACGC-CT-CTTCTCCA-T-TGACTAAAGGTCCCCCACCGCAACTCAATAACGTACACCTAT-CCTATAT-CCGT-A-GAACAAACAGATCGTCCCTCTTTCATTACCAACCCCAATTGAACTAACAACATATAGCCCT-AG-GC-CCACTATATCGCTGCAATACGTCCCTAAGCAGGCTAGTTCCAATAATTCAATCTAACTTAATGGCGGTAACACA-AAGCACAATTCGAAAGTCAACAGAGCCCCGCAAACACTAATAATAGAACAGTATTTA-AT-CAC-ATA-TCAAGACC-CCCCTAGCTCACTGCCGCACTAA-AA-TAAC-GGATGATAATAATTGCCCCAGAGGCGCTTACTCCTAAACTAAGCCGTAAGAAAAAGTACTTAACTCCTGTAAGTCAAATCACGGCCGCAAATA-C-ATTCCCTCCCCCCTAACCGGTCTTGTCGCCT-TTATAACGTAA--CTGCAAAAATCATCGCACC-GGCTCCCTGTCTGGTCTCTTGAGCACGTTACC-TCACCAGGACCTGATTT-CAATTCACTGACTGCCTATCCGCCTCTTCTCCAAGTAGTATCTCCTAAAAAAATAGAGACGGTCCCCAC-CGCACTA-CCGTCCACTGACCTCACTGCCTAAATCCTAAGCACCAAAAAACTCGTTT-ACCAAGCTATGCTAGCTACCCTTTGTTGTACGAAAACC-TA-TCTGTCGATGTCTA--AACGACTGATTAAAAACAACATAAACACTCATCACTTTCGTAG-ATACAAAACGGATC-TACCTAATCAAAGA-TATACGGCATATGCGGGT-ATGCGCT-TGT-CATGTAACAAC-T--G-CCAGCGCTTCAAGCGTCTGCG-CACAGAGGGCTTCTCATAATTACGAAGCAAGTTACCTATGACCCTACTAGGTAATCGAAGGGACC-ACATTCTTCC-CCA-CCGCTA-AAGAAAACATGTAATAACCGAATCATCAAACACAAACAGACCTACCAATC-CGTAACCCTCAAGACTACATAACCAACCGCTTAAAACAAATTAGTATATAAGCTAGAAAGGTGA-CCTGCATAAACTCGT-CATCT-CAGCAAT-C-CCGAATAAAAACCGTC--T--ACCATTGTACAAGCTAGCAGTCAGAT-TCGAGACAAGCATCCTCAAACTGTATCATAATTACAGATTAGCGCACAAGAGCACAAACGAAACTCGTCTCCTAAGCGGTGCC-T-TAAAGAATCAGCTAAATGACCAGCTACCAAACCCACACCTC-CTCCATGCAAATAAGTCATAAACTATAGCCTGTATCGTGATG-AACG-TACTTACCGGTACAAACGATT-TCGTACATA-CGTAACAGTAAATTCCA-TTCCCTTA-A-TACCT-AAACCATT-C-T-CAATACAACCTCTTGCGGCA-AAC-CTAAAGTGTC-CACC-TC-CATTCG-CCTTG-TCTCACCAACCTTAATGGTTCCAATCAACAAT-CCTAACCTCCACTAAAGTAAATCTCAATCTTCCAGATCTACGCGATATGATAGGCAACGAATCTTCACACATTGGACATTCGCGCACTCTTATATTGCCTAC-CGACGTATTTAAAGTCCATGGATATCCACTTAAAACACCAATCCCTGCTCCCTAATCAGGCGCTAGAAAATAACCTCATAGAC--CTA--GACACTCC-TATCACAAATAC-GCATGAACGTTAT-TACCTGACCG-CAT-TCCTGA-CCGGAAATCCATGATTTATACGTTCAAATCGTTCT-CAC-TTGCC-CAGTTACACTCCGTACTAGC-TACGAGCCAAGA-ACTAAAACTCGGC-AACCGA-CTCCTATAACTCCCTCTTCTTC-TAGAGAACTCAGGTACGCATGGCTAAAAAACA-ATGCAAT-ATCATTAAAAATA-CTCTTTAAGC-CTGACAAAGAGAGACACAATCCAG-G-TGAC-TCGTTACT-C-CATGTGA-GAACTATTAACACTTACCAAGGCT-TAAAGTACTATCCACCTATTCCT-CAGGAGCCG--ACC-TAAGGTCCAC-CT--GTATTAAATCAC-CCATTTGGGAACCCACCGGGAACTCGCAGTTTC-AACATACCGCAAACAAAGATT-CATAAAACA-ACCAATAACCT-CCAGGAATCAATCCAATCAC-TGTACAAACAGTCCACACGAGGCGAGCAACTATTG-AT-GG-TATGACACTACATTTCGAGCATCGGGGCAACAGGCGAAAGTCTCATCATAAGAACACTTACTTTGCCACCACACATTGGAAGTAATTCTAAGACTCCTAGTAAAGGTAGGATCTGTTGTCAACCATTGACGTGACGTGAGTCAAAGTTCAA-CAACACCCACCCGTAACACAATCTCCTGTAAGCCTACCTGTCAGCTGGCAATCGATTCATCCCTAAATAGAAACGCTAGGAACTTCTCATCAGCTATAACACTTTCGAATGCACA-CCGGGTGACCTATATTGGTTAAGTACTCATTTCGCCACAAC-GCCTATAATCCGTAGTAATCTACGTAACTAGTCATCCATATACCGATACTGTTGTAGACA-CCCGACTCTCCCTCTCTCAG-CCCT--CCAGTC-ATGAAACTAAACTTTACGAGAAGCCGACTCAATTGTTCTTCT-CCTACTCATGTATAATATTTCTAGATCACTCAGTCCCTTACAAGGTCCCCTTAAG-AACAGAGCCATTCACACCTCAGCGCAAATCTTATGAGCCGATGGCATTTCACGTGCAGAGAACCCGTTACTGCAAATCCCGCCCAAC-TCCGCGGTAAGCATAAACCTCCAATCGTCCTCTCACAGCCCCTCCAAATCACTCCAAGCACGTTAAGCCAGTAACGACTTGTAACACTTCATCTCTTACATCGAACCATAAGAAATTCACATACCCTTGTCAATTACCGGAAAGACACGCACTTGATAGTTATATAGACTCTTAGACGAATACATTC-ATCACACATTT-TTATTCCCTCGCAGATCAACCACTCACCCACACTCCCG--ATTAAAAATTCCATACCATCC-AGAGACCGAT-CAACCTG-ACGTGATCCAATCAA-C-CA-CACCCGTCC-CT-GACATCCCTATGTAGC-ACCGACCGGTCAACA--ATA-GGATGTTGCTTATGGTTAAGGACACACTTGAATAT-CTGCCTCA-CAGAGGTATACTAC--TTGTC-CATCGGATAATAGATGTATACGAATAAG-AACTATCACCACTG--C-GACCAACACACAACACATCGATAGAAAATAAACACCAGACATACAAGC-GGCGACAAATTAAG-AATAGCCGACCCCGATAAAAACG-GCGAAACCAA-AATGCTTCTTCTCTGATTTATTTCGGCATCAC-ATGACC-CCTAAGCTCTAAACTCGCTATTAAAGTACCTATGAACTCCGCGAAAATGCTCCCGCCCCGCCCAAACAAC-CTCTTTCC-T-GCAC---CGGAG-CC---CCCCAATCTAGAAGTTACCGCTGGGTAGCGACTTTAGGGAACCCACCCACTTACCACACACATGCC-ATTCGTTGCCCACCGTCGTGTTTCCCAACCACCAATCT-TCGTAAAGTGGAATATCGAAAGCGGCACTCTATTATTCCTCTGTT--TCCGCATAAGCAG-GA-GGAAACCACGATAAACATAATAAATT-TGTTGCACACCACCCACACAATACCCTTCTCACCCGACATTATATGCCAGTAACACC--GCCCTAACAGACCCGGGCTAGAAAATACTTTC-GTAGAAGCA-GACGGCAGGACCTACACGGC-C--ATCGCACCCTGA-AACCCCCCAGC--TGCTTTAGTACCTAT-ATTCTTTGCGTCA-ACAACCA-AACTATACGT-ACCTTCATTGCAGAAT-TTCC-G-ACAGTGCAAC-CTGAGCAGTAAGAATCTTGTAACTCTACTGAAAACCCAAACGGCCATA-AG---CC-GCTTC----CCACAAGATTTTTTACTTCCAAAAT--TC-GTCATC-CCCTTACAAATTCCGCCGCA-CCTAACAC-AC-ATAATACAT-CAAATCGAACCGTTAGCGACAAAAGTTG-ATCGACAGTAGAAATAAA-AAG-TACT-TGCTCAACATAGATTC-C-AC---CCAGATAGCATCACAATCCAGTTCTG-AAAAACC-TTCACACCTA-C-AGCTAAAAACTGGTGGTAGCTGTTGGATAACACCAAAAAATTTGACAG-GT--TTACCCGCGTAATACAGAC-AAGAC-ATTTGGTAAACAGACTGACTAATGCTTTTC-CAAACA-TAACA--AACTAAGCTTGATGTAACATCTTTCT-CAACCTATCCTTTCTAGCTTCTTCCAAGTAACAAATTTTGACCGATGCTATCTCCTCTTAAC-CAGGACCGAAAC-TCTTAGGA-ATATAAGCTGAGGAAAATA-ACAATC-TTACAACGCATCTGCACCAAGGACCGAGCCCAGCGAGAAAATCGTCTTCAACACGCACCCCAATATCATAGACTTTATAGCAAGAGACCCTGCGCCCTCCCTTATTTATAACTTTAC-GATATCATAAAATCAAG-TGTTATAAAACAGATACC-ATTTTTACCTCAAT-G-TTTCAGACTCA-CTA-CACCACTAACGTTGCA-ATAGAAGATCTTGAAAAAACCTTGCTTCAACTGGATTATTATTCAAAATC-A--ACATA-A-GCCCA--TTGG-AAT-CG--CCTACGAAGACATATCT-C-TATCCA--C-CTCCGCCGAATGGAAGCCAACCTAAGTAT-CGCCACC-TCATACCCCTAAAAATGTTAAACATTCCACCGAGGTTGCCC-GCTCCATTCTAAGACTTTTAAGTT--TAGGCCCATG-GCAAGGTA-TGCTCGCCGCGCCAC--A-ACGAACTTGGACTATAACGGTCACCGACG-TATATCAAGTCAGACCCA-ACG----ACTTCGC-GC-CCCGTATTGACAACGGAGGCAAATTCCACGA-AGCTCAGCTGAGCACAC-AACAGTAAACGCATA-CGCTTTATTTCAC-CTACAAATTCCTGGT-CAC-TC-AGG--ATATCTAACCACATTAATGCCAAAAATACAACCAT-ACTCC---GACA-CGTGAAACTATCCACAACCCCCTCCTTTCACTAACCG-CCCTTTAATACCGCTTTAAAGGCTACAATTAGACAATCTACAGCATCATTCAAATGAGCATTCTC-CCTGTAATACTTTCGA-GCTCAACTACCC-----CGTC--G--ATCCGG--AGTGGC-AAATC--CTTCC-TTGCA-CTA-ATACCATTCAAATCAGATAACA-CAACTCCTATAGAACAACACTA-AG-ACCTTTAACTTAAGCACTCATAACCTTAACCTTCACTGCTCCAAACTTCAAATCTCTTCCAGGACAGGAACTATTGTAAACGAGAGAGTTTAAAAATC-CTATTATCTTAAAGCACGACATTATCCTGCGCACTGGAAGTATCCAACTTGGATCC-CAA-AT-TAAG-A-TTCCATTCA-TGA-AGAATACTATAA-ACCAACC-ACG-AAGAACAATA-TAGGCA-T-CTAAAAC--T-CCCGGTCAACATAGCC--GAAACGTCAAGGCCTAACAGACCTTCTCCGG-AATTCTCTACTATACCAAAATCCCCTGATCTAAAAAATTAAACTGAGCTTAACA-AAATGGCAATGT-CTTAGACTTACC-AAC-ATCATACAACTTGACCAGAATACTCTTTTCCACAAAACCGATGATAATCATCGATCAGTCCACCCGTTAGACAAACTG-AAACTCAATAAG-A--AGCAGACC-GC--ACCCTACTAATGGCTAAATTATATAC-TAAT--ACTTGATAA-CTCA-TCTTTTAACATCGAATCGAC-CCGTTGCAGC-CCCA-AGAATATCACTATCGCCAGTAGCTTCGCCAATTTCTTAACGTTG-TAAAGT-CCGAAATTAGTATTTAATGG-CATC-AAGGTTCGAAAC-ACACTCGCTATCCGC-GACAAAGTAAAC-TACTCACTAAATCCT-CCAATTATGTCTTTCGACAAACCAGAATTTACG-GTGAATCACT-CAAGCCCTCCTCACCAACGC-TAGATAAAATCGGC-ATCTAGACGGA-GT-CTCCGTTTAATAAAATAACGGAATTCGACGTTCC-CCAC-ACATACCAA-GGGCAACGCTTTACTCCCAAACATTGAAGTAAAAATGGGGCTTCTTAGCATC-CATAGCCAGGGCGAGTCCCCTAATACTCCTAGTGAGAAAGCGAAAGGTAAAGCGCGTCAATCTAGCGAGAATACCCAAACTGCTTAGTCACGATCTTTGTCCTA-CATAATTGACTCCTATGGCATCTTAGGAAAGCATCCGACAGAGGATCTATCCTGCTGATGCATTGCGTTCGTATTAACCAACAATTTACCTGTGAGC-TCTAATCT-TATCATATGCC-CC-AACA-AAACTTAATC-GGTTGCTACATCCA-GCAATAATTTTCACCTCACATTCA-ACTGTCCACGTATGCATTAGCTCTATCCCATTAC-CCT-AAAACCTATCGCCAAT-CTT-AAG-CGTCAG-ACTCACGTAACCTCAAAACTGACTACCTCATTCCAAATCCAAGACTACGATCCTCCCACACATGTTAG-CCACGACGATAATAAGCGAACCAACCGTGGAAAATCACATTCGACTTTGTGC-TC-ACCTCTACTG-ACCAT-CGACCACGTACACTTCTGGTCGTCAACTATGC-GGGCCCTCCGCCAGACTTACAAG-CAACTT----CTTCGCAGAAACTGACATTAAGAAAAATTAAATACAGGCGCCATACGCCCCAAGTTATCCATTCTTAGATTCACTTTATT-CCATCAAACTATGTATTCGCAACTCACAAGGCTC-CCAATGCTGAC----TGACTACTACCCAGACCTTCATAT-AAACTA--CGTGCAA-CTG-TAAAATGCTTCACGACCATTCGCTTTCCCCACAGTCCAGAACCTGA-CC-CTA---AGAAGGTAATATAACGCCAAACCAT-CGCTTCAACAAAAACATAATAG-CAATTCCAAAATGA-ACA-GGATAGATAAGTTAAAGGTTTACAGTCTTCCCC-ATAACTATGTTC-CTGATCC-TTAAATATAATGGGTTA-GTA-C--AAC-GCAAACGCT--CAAATCGGATCATTCTTGGA-TAATA-ACAGTCCC-ATATAC-T-TTTCA-GAAATACACACGTGAGGGAAAAC-TA-GAATGCCGT-AATATAGACACTTTC-T-GTCGAG-AAGAGGTG-AT-TTCG-TTCATCCA-TCCTTAAAGTTACCGCAAC-CGCTTCTCTAGCC-GATTC-ACA-G-AC-C-CTTCCCCCCAAGCTTACT-AC--TCC--CTGATCACTACCATTAGCTAG-TGTACTGTCAG-GCCGGATTCATATAGTCACCACCACGCTATGGTCA-CCCCTAACACTAACCTCCTTAAAATACCAGCAATAGGCCATACTCCACGGATAGATTACCTAAAAATCAC--ACCTAACGCGAACAAACTAGC-CTATTATTCTC-ACAATTGATTATTACACCAGACATCCTAAATTTCTTCTTGCTAATCCCTATCTCTCGACGCTAGCAAAGAGACTACTCC-TAACGTTGTTTCTAAACACCTCGCCAATCTCT-TAAA-CCGAAGTTTGACAATACCGTTTCTCAGAACCTCTTCTAAATATGGCAAGCTAAC-TTTCATCTTAGAGTTATCCAGG-TAGTGAAAAGTATACGTTACGTGCTATATCAGTATCCGCAACTATATACTCAAGGATAATCTTCC--AGCACACACTTA-TA-TA-CACCTTGATTTATG-TA-TCATAAGAAATACGCCACAGTCTTGAA--C-TGAATCGTGTGAAAATTCTACT-CCATGAATCCAGTGACCTTACGAATATCGCGAAA-TCATCA-CTGAGCGTTCT-TAGACCCTTCTGGCACCTTTATGTTTTTAAAATACT-AAACATGCAA-ACA-ACTTC-TTATTGTAAAGCCTAGCACTTCAGC-CACGGATCCTCCG-A-AAACACC-AACT-AGATC-TGA-G-T-GACAACATGAACTCC-CCTTA-ATAGCGATTGAGCACTCTCATT-TA-CC-G-A---CG-CAACTTTTGGATCAAT-TTATCC-CCCCCT-GATTCAATGACC-TCCGCATCTAGCCACCGCT-T-CTA-A-CA-TGTGTTAGTGTACAATAAA-AGACCCGTACAACCCAAAAAATCTGCTCA-CTCCGCTATCAGCCCAAGGACATCAGCGTCATGG-AACCCCATA-TCCAA-ACCTCAA-TATAACACA-CCCCTCGTCTCCGTAATAATTATAG-TTACTCGGATTGTAAACATCCCTCGCGTATTCTTCTCAGTACAGTAATAGCTTATCTTCCCGATGTGCACCCTCCTGCGCGGTTTACCTA-AACCATACTCATCTTCCACCAATCTGGTACATCAATTCCATAACACCCGATG--AAGTTTTAACCGCCTGAGTGTTATTCAGCTTCAGCATGTCAAACT-TACCAAATGCC-GATTGTG-GAGCTCACACATTATTGAACACAACTC--T-TACTG-ACTTAATC-A-TTCTACCCTACCGCCATAAAGCCT-GTTA--CGAATAGCGCATTTCTAACATAGAATCACCAACCATCTAACGATGATTCTCACCCCATACCCCATCATCTAAGTGCCACCCTCTAGCTAGAACAGTTCTTCCACCGCAACGTTCCGGGTCGTACACTTAACGTTACAAACGAGCCTAAA-ATTC-AGGGTATACTCTGCCCTTT-GAATAA-CTAGTCATA-TTCAGCTCTATAAACCAAAACATTATTTTAACTTTCCCTAAATTTACGTCAAATAAACATT-CCTTAT-CGCTTACCTCTTTAGTTACTGCTCGACGT--A-CGCATGCGACT-ACA-CTACGGGGTACCCAGAAAAAAGCTCAACCACGCAATAACA-AACTGTGCATGAACCGTTAAACTTACTAGAAATAGGTCTCATCCCTGCCACTAGGATCCAGA-CAAT--CAGA-AT-AGGACACT-T-CCCCTAAATCAATT-GTAATCTCCA-G-AGGACTAGC-CAAAT-C-ATGCCTCGTAA-AACTAAATCATAGGCAC-CTCCCTCGCACC-GATACTCTCCCTAAAAAACAC-TATC-CAGATCCCCAAC-TACATAA-AACTAAAATCAGAC-GCCTAAAACTAATATCTTCTCA-ACTATAAAACCGATGAACACAC-ACCAAT--CACTA-C-TTT-CCTTG-A-CATAGCCCTC--ACCAGGAATTTTTAC-GCAAAT-GAA-TATCT-CTG-TACC-ACTGACTTAGAACGAT-TCGA-C-ACAAGAACACAAACAAATTC-ATT-TAAGTTTTTTCGCTTTAC-TTTCGTATACTGTGAAGAC-TCCAGCCCTGTAGACCTTTCCCAACCATTAAATCTAGTC-TTACTAG--A-CAGTATAAGACTATGGCA-TACGG-TCCTATAAATACTAATTTC-TGAACTCTATCAA-TTATACAACTC---ATATT-ATGAATGT-CCTGAC-CAGTAC-CGA-GCCAGGCG-AGCGATGATCAGAGCGTTAAATATCTCAAAGCCATCCCAATTAACGCCCAAGAACCCTCACTAGAAATAGTGCCCCGAG-TTCTCCCCACTATA-AC-AT--CACCAC-T-ATA-CA-CCTAGGGATCAT-ACATATAGACAG-CGCTCAATGGTACTACTAAAAGTCGAAACGTTT-TGAGAGA-ACAAA--ACCA--CATGTAATATCCCGTGGATTC-TAAA-TTTCA-GGCCAACAG-TGT-ACACAA-CACGGAAACTACCACTAATATCTTAATAAAAACCCCAGAATTCAGGTCTAAGAAAACATAGATACCTCATCCTCTCCCTCCGTGCT-TAAG-GCAATTA-GCATCCT-CTGACATATAGAATTTAATATTACT-AAGAAAACGACAACCCAGGCGATTCTATTGTACACTCGACTATACAAAGTGTTTAGC-TGACTCTCGCGCACCCGTATCATTCTCTCATTCCCCTTTCTCAAAAACACGTGCCCTTAACCAGCCATACGACAAATTTTCCAACATTTAATCACAAAATGATCTAAAGACAACCATGCATT-CGAAGGATAATCCCGCGCCCTCATTAACCT-TAACGCTCTCTTACTGATCCCACGACA-ACGACACTCTTA-C---TAAACTCAACTACAG-CCTCTACGCTCCAAGTTAATCCTTACAAGT-T--TTACAAAGTAC-CGAAGG--TACA-ATGACCCAGC-CTTTGGTATTCCACTC-AGACTCCTCCAAATACATATT-TGGA--ACAAGCATCAGAA-ACCATTCCCTA-AATCTTTCTACAGGGAAAGAACCCTAACAGG-GAATCAAACACC-TCGAACTTTGTCCTCT-CTACCGACGGCACATCGAGCATCCCCGATACCATACAGCACT-TGCCGCTCATACCAAAGAATCGCCGA-ACAATCGGGTCCCTAATTACT--ATTTAG-ATCCTCCAACTCTCTCTAACTC--C-CAATAGTAGCC-TTAGTGTAACGACTCAAATAACCATTATATCAAT-A--AAACAGTGCCATAGAAAAC-TATTATACCGCGTTTCAAAACCTAAAACCCTGGA-CATCCGCGCCC-ACAAGAAATTCCAGTCC-A-TCCACCTAACAAGTCTA-CAGA-TAATAGCTAC-ACATAACAGGCACGACATCAAAAACAAC-CATTATC-CGTGCATCTC-AGATA-TTCCAGCTATATCAAAAAATGAATTCCCGTTCGCACTAACTAAAACCACCTGTTTAATTATGTACGGGCTCCATTCACTTAACAACAGAGAC-AGAACTCT-CCAGAAGGTCCT-TGCTGGAATAACATTTCAAGTCCC-TAACC--ATCCCAACGAAGATTCATCAGCAC-TCTGG-CTGAGAATCTTCGAGACAGCATT-TAATCGGCCGTCATTGC-ACATGCGTTCTATCCATTTGACACAAT-CACACAACGACAGCCCT-G-AACTTCCAATTCCATAAACCCTCGTGTGC-GT-TAC--CTTGTGATACGAGAAAAGTGAAGTCATCG-TCA--AG--CCGCCCT-CGATCCGCCAACATCGGT-ACCGAAG-C-ACACCCACACCTCACGGGGTC----GCG-A-T-TCTAACCCATCTCATCCAGCTA-ACTCTAAAAAGTTTAACAAA-TT-ATAT-AG-ACTTCC-CGCTAGTAACAAGAAAC-CGTTGGAA-T--TCATAGATTTAATATCCATGC-CGA-TTGTAAATGACT-ACTC-ATAACCTCAAACAAACGTTTTTTAGCCCTCGCTAGACC-GTTAAAACG-CCTATGTCTACTA-AATCCAAACATTAA-CTTAGACTACGCTAAAATTGCATTCCCACACAGTCACCTCCTGA-G-A--AGGATTGCTAAGTTACGG-AGAGTACG-T-GTCCCGTGT-AT-GAC-TCCGCAACAGAGC-C-TA-GTCACCC-TG---CTC-CCAAATCACCCAAAATGAGACA-AAGCATAATACAAATTGCTTTCAGATT-GTCAATACCGC-ACCGAAGAATGCTACTCATAACCTATCTAAACTATTGAAACTTGGTGTTTA-TGAAATATGGCTCAGCGACCCCTAATGCTAAGGCCAAAATGTCT-AAGACTTGTCAATCAAACGGCATGAAACCAAT-CCC-CT--CT--TCAGA-ATGA-CC----CAGTAGTTTAGG---A-CAG-ATCC-TTGCC-CTTTC--ACA-TCATACCCTAC-CAACTAGAAACAGGAGGCAGAGTTGTATAGCAAACACATAAATTGTCACCCAACCATCCTTCTAACC-CACTAGAGGATACT-GTG-CGGATTCGCAAATCCTCCTCCCCCGGTCTGTTAAACCCCTAACTT-AGAA-CTCTTCCATTACTCA-GCTCCTCAAATACGACC-ACTAAACCCCTTATCACGTCTGC-AACTTACCGCAATGCC-GAC-C--ATCCGGCATATCAAATTCAATAA-AG-ACTTG-CCCACACCAAACGCGTTAACCTAACACCCTCAACACCTTGCTAATTTCATTAACGACACAACAGAAGCCTCC-ATAGC-CG--CT-AA-CGA-CCCCCGGGTAAAGTA--GCTA-T---CGA-CG-TG-G-GAATAAACACTTGTACCTCAAAAAAGTCCAAATGCCTAACCATTTTCAG-TACTTATACTGACAACAAATAGTACTCCCCCAACAC-CCCTTAAACTAG-TTACTAACTTCACGCATTATAATCAGAAATTCTGGTCAACAGAT-CTCCATATCCTAATCCGTATGTAAAATGCATAC-CGAACCATGTTCTTCAGGGCACAGTAACA-TCAGAATAACCCTATATAACTAACTTTCCGTAAGTATGTGTGGTCATCTAAAACCATAACTTTCAATGCGTGTTACTAAAACACTTTAAATACGCA-ACGGATCCAGGTTACAAAAACACGGAAAGAACGATTCACTACAC-AATCACGACC---CCTC-CATGATCAT-TTTCCCACCCGCACATTCAC-G-C-ACCCATATT-GCAA--TAGTACC-G--TTA-TCAATGAT-TTCTCAAAACGACCTTCTTGATAGTCAAATTACAGGTCAACCTATTAAATAAACA-AATACAA-TCCAAACACGCTTCACTCGCTT-CAACCCAGCAAAATTACGAACCCCGGTCCAAGAATCTCATTATAAATTTCCACCTCTCCATGTTCAACTATTCTAAGT-TTCAATGAT-GCCCTACACCT-AACTTCCAAGACTA-G-AAACCCCAAACGAACC---C-CCTGATCCCATTTACCTCAGAATACAATTCGTACTCCGAAG-C-CGACAAAGACA-GCCGT-ATAGC--CCAGCTCGCAGAAATGACATCGGCGAGGGTAAATTGACAAGGGTTC-TTA-CATTTCAACTTGACCCC-ACACCCGCCAATTC-CTCGC-AATTCGCGCACTCATATTCGTACACGTA-CCAACAGGGTTATATCCATTTGTTGAAAAC--CAGCTAC-TAAAAATAAAC-CCCTA-CAC-GCAGAGACTTC-TATCGCC--CACCACCTACATAACTGGAATAACAGGC-AAATGCAAAACTCTCC-GACTTTCGGCATTCGC-TGAC-C-CCCCAAATTACCCATAGAATTCTCACTCTGACGAACACGCCGCGCCACGTCC-GCACATAACTTAGACGCAATAGGCAGACTCTAG-CTCTTACATACAATAACCATTGTAATTAGATATGGCTCAACTG-TTCTCGAAGGACA-AATATCGCCGTTACCCTATCGCA-TA-CTCTACTC-GATACCACTGGATTAT-GA-AAA-GTGCGCAACGCAATATGAATTAAACCGACGTATTCTCAAACTCACCCAGCC-CA-CAC-TAATAGACAAGC-GCAAATT-T-AACCTCTTATTCTC-CCTTAAAGCTGGTTAATCCCAGCAATAGACACTCTAGTAAACTGAGAACGACCAAGACCCGTTCTGGAACAGAGATCCACCTTTCATACCGCCATAT-TCAACGCCCGATCACAATC-GAGTAGTTTGCACTCGACGAGGTAGACTATCTTC-CA---CCTCCTTGTTTAAC-ACCTATTATCGATGCCGTACCTAAAGGTTCCATAAACTTCAGAAATCCAACAATGAC-TCTGTATTCCAAGCGCTGAAAACAAGATAGAGAATCATAC-CCAGTACATCACTCTAATACAGAGCATACGTCGAAC-CTAT-AGCCAGAACTTCGGAAACACAC-TCCTTCGAACCGTTTACCAACGGACGCGTTAACT-AACTTCTT-TCTAAAACCTAC-CAGCCTCTCCTATCCACTAATA-ACGATGA-GTAT-AGATAACTAAC-GCGAAGATGGCAAACACACAG-CCT-GCGGCATTTAACTTGGCA-CACGACACCT-GCAATCTATATATC-AAGCTCCAACGTCAGTCACTGTC-TCCGCTCTCAGTACTAATCC-CACACA-AACT-ATCAC-TTACGAAATAATCTT-TGTCAAATACGTTGTACTGACTGGTCTACCTTAAAGTCTCTAACCTAAGTTTATCCCATCCTTGGCAAACTCGTTTTATACAATTAACCCCCACGAAGATGTGATTACGTGGCTTCGTATAATCTCACCGCAAATTTTCAA-AAAAATCCTGTCCTAGGAAACTTATA-TCGTGAACATGCCCGAACTC-TCTTCGAGAACCCCTATTACCAAAATTAAACACTAACAAAAGTTCATATGATTGTATTCAGACATCCCTTCTAAACCAAGTAAAAAACCCACTTCCATCAGGGAAACACATTTAAACTGTGACCT-ACAAACCTTCATACTTGTCGCTATACAAGAGACTAGCCATACACCAAACATCGGTCGAC-AC-CAACCAATCCCTACATGCCTTTTAACCCACCTAAACTGACGACAACCGCACACGTCGTCGACACTTTAA-CTTGTCAATCGCAATAAATACCGAGCTGC-ACA-CCGA-TATCA-A-A-GAAAC-AGAAGAG-CGACTCGG-CG-TCGCT-GTCCCCCTAAGATT-A-AA-CCCTTTCTTCTTTAGCTTAGCTTCTACGACCCCATAATGAATTGTTAACCTACTTTTCATTACAGTTAACCGAC-AAATACTAAAGCGTAGAC-G-ACCCGCGCAAGTCACTCCCT-ACACA-ACTCTAATACATAACCT-TTTAAGCTACCGTCTCACAAGTTGCAG-TACCACCTGGATATCTAGCCGGTCTCCTC-CGCACACAGCTCGCTCTAACTCATAGTGAATCAAAG-TACCGTTG-CCT--CG-CCGAACTC-CACAACGGACGTTAAGTTTAGCTC-C-CC-ACAC-GCTGACCCCACCCG-CATGGCGCCAAGCTTCAGCAAAAACGC--ATT-CC-AC-CCC-CATGT---TCAAAATAGT-TAACCA--TTCCT-AAGCATAAACATTCTAATGTGCCAAATGCATATTTCAGAAAAC-GAC-CCCCAGAATTAGATACAACAACATCACCTTTGAGCTTCCCATAAAACTAGTCTCGTCACTTAAGGAATAGGATCTTTAACAGATCTCCCACGTCCCCAACCACCCAATAATCTGCTCC-AAAGAGCACCCGACC-T--ACCTTA-AAAATTCTGGTGACATTATAC-GTCTAAACGAATGGCCG-GACACACG-GTGACAGTATGCAGAATTAAATCCACACACTACC-----CTTCTCAAGCGTA-TAATAAAGTAC---CG-TTTACACCCCCTTCCACAGAGATGGAGCTCTCCGAGCCATTTATAATTGCA-GTAGC--TTCATCT-TTCCTA-AT-GCT-ACACAG-TA--TACG-T---CGACGATT--ACA-CCCA-A-TAAGTCAATGGAACCT-GCCCTCA--ACGC-CTGCC-TA-T-ATTC--ACTC--CCTCTAGTC--CAAGAT-GCACACCTAAAAGACCACCCCAATACAGTCACACTACTTA--ATA-CCACATACCCACCCATG-AT-GTA-CGTTGC-CCCAATTACC-ACACTCTTGACAAGTGGTGACAAGGA-AGTATGTAGGCCATTAAATTTAAACTATGGTCACTCTACAGTACAATAGCATCCTGTTACCCGAGTGGCATATTACTAATCGTGTAGCCTACCCTGCCACTCTCGA-TCTGTATCTAACAAACATCCATATCAC-TTCCGTATACACGG-TATCGTCCGCGC-C-AAACT--CCA--A-C-T-CATAAATTATGC-CTAGA-A-CACAATATTAC-T-CTCTTTTTGGCAAATCGCCTCT-TA-ATATCACCATCT-CATCTTATAAGTCATTCATCAGATTAT-AGCTA-CGTACTC-CTCGT-ATATCTTTCACTCCGCGACCCCAATGTAGTATAC-TAGCC-TACAACAAACTACCTCTCTATCATTATCGCGAGGAATACTG-TGACA-G-TT--A-TAAAAACACATACTGAAACTT-GGTCGTACCTTAC-CTTGGGCCTCAGACAAACAACAG-GAAC-TGAA-GTTGATCCAC-A--GATTACCGGTCTC'
    '|   |:|||||| :|::|||: |||   ::| || | | :|::| ||  ::|| | || :: | |||| | |||| |:::|| :  | : :::|: || :|| |:|| | :| ||| : ::| || |:| || |||| |||:| ||: :||:|||  | | || | :: |::|||  :  :||: ||:|| ||| |||::| |  | ||:| |||: ||  || ||:::|  | |:|  ::|::  |:  | |  ||||  |:||||  |::||:: |:| :|| |:|:||| |:| :: :|||| | : |:|  | : ||| :|:  |   || ||:|:||| : |||:: :|||| || |  || : |||  |:|:||  |||:||:|| ::||:| ||  : ::| |||||||  :||:| :| :| || |  | |:|| :||| :| |||   | ||||| |::| ||:  :|  :  :|:::  |: :|: || | || |  ::|| ||| | :|| :||: |:| ||| |  : : | |||| | | |||||| | :|||||:| : ||:|| :| | |: || : :|:|||   ||:||:  ::|: || | || | ::||||||  :| |:|||  |||: |||:||:   ||:||  || | :  || :|  |||  :| ||  |||| |||: ::||| :| |: :||::|::  |||  |:| || ||:|:|: :| ||  :: | |: :|: ::|:| | :| ||| | || ||:  :|:| :| | ||  | | || || | ::|   | ||||| || |||::|  | ||| |: | ||:||| ||:::| ||||   | ||  || | || ||  | || ||:| |: || || || |:| || ||: ||  || :|||| |:| :| || :|||||||  |||| | :| :|| ||:::|  |::    || |:||||| | ||| |  ||||| |   |   :| ||||| :|:::: ||  ::||:| | || ||  |: : |  | |||:|| ||||   | |:|::| :||:|||| |||  ||   ||:| |||  :||: ||:| ::| | ||::: |  ::|:| |:|   || :|| | || ||| :| ||||::  :||| | | |:| |  |  || |||| : ||:|| ||   |||||:|| |::|: :| ::|:|| :||: | || || |  ||| | :||| |:||:|  |: | :|:|||: |:  :|||| |::||||  |:| ||||| | ||||  | :| ||:| |||:|:||: ||  ||: | |::| ||:|    ||| : || |||:|||:: :||::: |:::: |||| :  |:|::::| ||| : |||: |    |:|:|||| |||| |||:|: |:||||:| |: | | | | | :||:|| |  ||: ||||  :  | :|||| |: |:|| |||  | |   |:|    || |::| || :| ||  || |:   | |  ||:|:||||  ||||   :| |||: :|:: || :| : |||:| | |:|   :|  |   :| | ||| ::: :||| | |   |: ||  ||| :::::||||| |||| ::|||  :|||  |  ||: :|| :|  |||| |||  | :|||:|: |: |  ::|| | |||   :|:|| | :|:: | || : :| | |: || |:||: |: | || | ||||||  | :||| || :  :|| ||| |||  ||  |  |||:|:| || |: ||  |||:  |||||:|  :  |  |||  |:|| : |:||| :| ||  | : ||||  |:| ||  |:|:||  ||| ||:|||| ||:: ||| |::|:|: |||||| |:|| ||| | | | ||||:|||||:|:| |:|:   |  ||||:  | |::| || | | ||:| :||::| |  | |||: | :||:: |:||| || | |  |::||||  ||:||| :||  :|  | | || ||::|||:| | || |:||  : | ||  ||||:| :||   | | :|:| :: :|::|:| | :||:||   |||| : ||||::: | | ||: ||:|| |:| ||| :|| | | || |  ||:|| || | :| |||: | |:: : || || ||:|  ::| ||| ||| | |: || |::||  :||:||   ||  | |  :||||||  |  | |:: |||||| | |: : :|:|:| || :|| :|||  | || |: :| || :||  |||| |: :| | | ||: ||||  ||| ::| |:| :|| | ||||| || |  ||:|||:|  ||| || || ::|| ||:||   | :  :| || |:|:| |:|:|  ||  :: | | | |||||: || ||: :||| || | | |: | ||:: ||:: :| |||  | ||:|:|:| |:| :||: ||| |:|| : ||||  :|| |:|| | | | |||||:  | |: |:||:  :  |:||:| ::| || | :   || || || :| ||: : |:| ||:|||:| |||  || ||::| | |||  || |||:|:| |: || |||  :| |  :|| || |||| |||::: ||||  ||| ||| : ||::: |:|| : :| | |||| |   ||: | ||:|  | : |:| :| | :|| |  |||:| ||||||| :|  | | ||  |: :|:| |:| ||||| :||  ||::| ||   ||  ||| |:|| ||: ||||| || :|||:|:| ::|| |:||  |||: ||| :  ::|  | :| | || | |  || ||:  | |:  | :|| | | :| || || :|| |||  ::|||: | |:||:| |||:| || : |:| | :| |||| | |:  ||||: | ||  :| | | |: || :|: || ||| ||:: | |: |||||  :||  | || |::||||||| | | : : || |:  |||   | || :|  ||||| | | || | |||:|| | ||:|||: |||||:: :: || |||| | :|| :|||:||  |: || | | :||| |:| | |||  ||  ||: ||  :| ||  | | | ||| ||  |: ::|: |||  | | :::||||| | :||| |:| | || | |:| : |||| |::|: :|: : |:|||: ||  |:|||| | :|| :::| |   || : ||:|:| || | | :||: |:: |||:| || :| ||| :| ||:| ::| :: | :| | || :|:  | :|:| |:|  |: :: ||||| | | ||||| |||||  ||: ||| :|  |:::| ||  :||    |::| |  | |:   |:| |:|:|: |  ||  |||||:||: |||| ||:|:  |:|||   | :|::|||| ||  |||  |  ||| | |  | |:|| |  | || |:::| |||:| | :| : ||| | || |  ||:::| :||| | || ||:||  |:| |:: |:: :|| :: :| |:|| : | |:||| ::| :|| : ||:|| || | | || | :|:|: |: :|||:||::|||   | | |: |::|| |||||| |: |||| | |||   ||  ||:||| | |:||||:::||||:| |::|| | :| |||| | | || |  |:| ||:| ::|||: ||:||  | |||| :||   ||||    |||||::|:| |||| | || ||||||  :|| |||  |||| : || ||:|| | | : || ||:|  |  :|||   | || ::| |||    | | |||:  ||:: | | :|||:  | :|| ||:| | || |  | ||||  ||  |||:: : |: | |||   || ||| ::||:|| : |||| |: :||:  :| || |||:|: ||||| |  ||:| |||  :| ||:::|:|| | :| :   :|  : |:|:||:|:|: | |::|: || | ||| ||||||| :| |  ||||:||  : |:  ||||| | :|: |   |  ||| ::    ||:||| | || |  :| ||:|  |: |:| | | :| || :|||| |  | |: |||: | |: :|||:|| |  ||| | ||| |||:|||| ||  :|| :||:   | ||||:| ||  | |:||||: |:||| || || || |||  |    |:||  || || |  |: :|::||:| || |: :| :||| ::| |||||||| |  |  | |: |||  || ||: |   : |:| : :||:||||||| :|: :| |||:| ::: | || |:||  | :| |  |||  ::|   |: |:| || :|: ||: |||:|  | | : ||||| |||:  |::| |  ||||| | ::|||:|:  | | :|  ||| | :|:| ||  |::: :||: |   || :| ||:| |:  ||| |||| | ||||   ||:|| ||    | | ||||  |:|  | | |||: :|:|:| : :  : :|||| |||| :|||| | |: |::|  ||||  : |||| | :|:: :|:|  || |:  ||:|| :|:  | |   |||||| :    :||| |||:||:|: |||:||:  :| | :|:||||: :| :|||:| ||:|||:||::    : |  ::|||| ||  | ||: ||  || ||: || :  |: ||:|  | |||  || ||||  || | :||  :||||:: :  |:::| ||:|  | :|||: | |:|  ||  | |:| ||:|| |  ||   : ||||: ::| | | :||  ||||||| |:| ||:  :||: | || || ||:: || :||:::|| | || |:|| :| |:||  | || | ||||| | :|  |:|||| :|:|| :|  |  :| || | :||| ||| ::| || :|| ||   || |||||     ||||:|:      :| :||:: ||  || ::|  | |||::|||: ||| | |  | |||:: || || ||  || :| |::|: : :||:    | |||:||:||| |||  | | :   | ||| | | | || |:|: |:|| |:||:| | |    ||:::|||| ||:|:|  |: |||:: |||:||| :| || ||:| | |||  :  |||  |||::||  :::|   ||: |:::|||| ||| || :|  :|: ||: :|: ||||: | :: :| |||||:   :|| :| ::||||:|||:| | |:::|| |:|||   :|  |||||  ||::||  |||:|  ||:| |:| |:   ::| | |||  |:||: |::| |  :| || |  |:|||| || :::| |:  |||  |:|  | :  || |::|::| |: | :  ||| || ||| :|||:| :|::|: |:||::::||| | |||  | |: |||:  ||  ||:| |:||  | |::|||:|  :|  |||: | | |  |||| |:||| | |: :||:::  ||   | ||  || ||:|:| |:| :| : || :||| ||:|| | : |  || |::| :  |||: ||||: |:| |||||||  |||:||| | |||| : || |:::||   |||||:|| |::|:  :: | ||:| | | |  :|:|: :  ||||  ||:| ||| ||  ||  |||::||| |:|| | | :|||     || | | |  ||||:|||:|| |: |||  |  ||| ||  ||||| |::|| :|:|::|:||| : |:| | : ||| ||  |:: ||:| || | || |||  || ||  |:: |||    | ||| || |: :||||  | |||::|:| | |: ||: | | :| |||| :||::||:::|: ||  | |||    | :||:| :|  || |||||  :|| | :  ||  || :| :  | ||:| ::| ||||| ::| :| || |||:|  : | ||: | || || |  || || ||| | | || : |  ||:|| :: |:|   |||| ||  |: :|| ||:|  ||||   :| | ||| :|| |:|  :|::|||||   : || |||: |: |||    | :| ||: |:||:: |||:||||:||| : :|| ||:  |::| :|| | ||:|| | ||::|::::|:::| | | :|| ||||||      |||  :  |:||:|  ||| || |::||  |  || ::| | ||| |:| |:::||::| | :|:|||  || | ||:|  |:||:| |:| || ||| || :| || | | ||:: |  |:||  |||: || :||:|| :| ::| | ||| |:| | |: |||:  :    | :|:  : ||  |:||||  |||::| ||||:: : | |::|: ||::  |:| : |||:|| |::|: |:| || |:|  | |||| | ||  ||| |  ||   |:| |   || ||||:|| | : :| ::||||: :|:||| | ||::| |  | ||| |:||:  |:|||  |||:|:||:| | || : | ||| ||:| |  |:::|  ||||: :|||: |  ||| ||| :||: : :| |:| | | ::  ||  ||||||||||| |:|| :| | || ::| ||| | | || |:| | ||| | | :|||||| ::||  || || :|| :: |: || :| |||| |::|| |:| |: |:|| ||: ::| |  |:|: ||| ||  :| |:|| | :  |:||   | : :| |::|  ||::|:|:| |||| :| ||  : :: ||:||::  || :::||:| | || | :|||| :||:|   ||::  |||  |||:| :||| | :||: |  |:|  | :  ||:| | : | :|| | :| :| :||| :|:| || || | |:|  || |||:|| ||::| ||||| | ||||| : ||||| ||:| :  | :|   ||:|:: ||||: :  :|:| || ||::  |:  ||: |||||| |:|:  :|  || | || |  :||:| || || |:|::||  ::|  :|||||  |:|  |:||  ||| | :|| ||:  ||| :| | ||| |  ||| |:|:  |:||::||  | | |:| :|| :|| |:|: ||:   | | : |||:|| ||: ||: ||  ::||  ||:  |:: || |: |||| :::|   :||:|||||||:: | :|: ||:|| :   :|||| |: |:||:|||  :|   || ||||: :|| | ||  :|||:: : |||: |  ||:  :||:: |:|:| | |  || | | | | ||| ||   | || :|||: |:|| |:| || || || |  | |:| :|| ||: : ||| | || |  |||:  ::|:    |||: |: ||| :|| |  || | |:| | | : ||||||||  |: ||:|  ||:| || |  ||| ||| |:|| | :|:|:|:||: | |   || | |:||| |||   ||| || ::|| :| | | :|  :| | |||| |  ||  :| ||:| :| | |:|| :| |:  |:||| |||:|||: | | : | |  | | ||:||| | |:| |:: |||   | | ||::::||: |||:|: | ||: | |  :||:|| : || :| ||:|||    | ||||| ||||   ||::| : : |||| :||:||  |: |||:| |   :|| :|:| :|: :::|| :|:|:: ||  ||| | ::| :||  |:|| ||| |||   ||:| | ||::||:||    ||| |:|:| | | | ||:|:|:| |::|||  ||: |:: |:| ||||:|: |:|:| |||||||||: |   | ||| | | | || :|  | |||   |:||:||: :||::|:||  :||:|    |||| |::   | ||: || ||||||| | |||| :|: ::||    :|  || | | :: |::| : |||| | :|||: :||   |:| || ||::||  || :| ::|  || |  ||| :|:: | ||  | || |:|  | |:|::::| |:||| || : ||| :|| :| | |::|: |:: | ||| ||| | ::| :|| || :|||| ||: :| :||| |: |:|| | ||  |: |:    |: || |:|| ||||:| | |||: |:  ||| |:|||   :||: : :: || ||::| | | | || | | : |||  |: || | | ||  :||  || :|||| |||: || ||||  :|:||:||:| :|  | | ||||| || : |:  :| ||: |::|| ||  |::|:|   ||   | :: ||| | |||| | : || | :| || ||| |::  ||:| :: ||||  :||||:|: : ||:|: |||  | || :: | | :|:| ::| : |:||| || |||:| ||   :: |:|||:||  ||  :|| : ||  ||||: ||:|   |||  :|| :|||:|: | :||::| ||||| |||| | || :||| ||   ::|:||| | |   | ||||  ||  || ||  ||   |: ::| : || |||||||::|    || ||:|| | : :||:|  :||  :: ||  | :|  |: :| || |:||  || | : || :||||||  ||  || || |:|:|| |: :| ||| |||:|  | : || ||| |: ||:|:|| |||  |||: ||  | |||||| ||:| : :: || |:  |:| :|| |||||:|| | ||| |: ||  :|: :|:||| || :||:  |:  :|| ||  || :|||||:|: | :     ||: || ||:| | ||| :|: ||: | | | ||| : |  ||||:||: ||| ||||: :  :| | | :::|:|| |: | :|||| | | : | | |||||::|| ||| || || ||: ||| ||| || | ||||| || || : |   || |||||:::| |:|| |  |:| |  |  |: ||:| :|| ||  |||| :||||  ||||: | | | | | || || |  |::: |||: |:| :||||   :|||   ||::|:| | | || || || :|  |:| |||| | | |||||  ||||  :|  || | ||||| | | ||| || : || | |||| |:| ||||| | |: ||:|| ||||  |||    | : |:|| | : ||    |||:||| :||  |:  :| |: ||| ||| | | | ||| | ::| | |    | || :||| |:| | | ||  :||::| ||: :|:||   |||:||: | | | ||  :|| ||:||  ||   |  | :|::||: |||| |::|:|| |||  : |||:|:|| |||::|: :|:|| || |:|| :|| :| |  |||  | |:|:| |||:|||| | |:| |||||: ||| |::: || | |  |  || :||||:|:   ||:||   |: || |||   :| ||:  || :|: : | ||| |  || |: || |||::|  | | | ||   |||| ||  |   |||| |:: ||   |:::| :: |||||||  :||: |||:|| :::|  |:| :: ||||: || :|| |::  ||||| |||| |||| ::|: |:|||  :|||   :|  || :|| ||   | | |: : | || | || : :|||: || :|: | |::|| :| ::||  :|| | ||: ||  |  |||| ||| | |||  :||:: | ||||||  :|:| |:|:||  |||   : || |||  | |: |  ||: | : ||||||| :|| | ::| ||  ||:||  | :||::| ||| |:||  |:|| :|  ::||| | |   ||::||| ::|| | :: | | | |  : ||| || |:| | | ||| |||  || || |:||||| :  ||  |:   | || || |:||||  | |||  :: ||| :||   |||     :|| :||||::  :|:|:| | |:|  |||:||: ||: ::||  | || || |||: :| | ||:|||||| |  ::|   ||:| | ||| ||| : | |::|:  |||  ||||| ||:: :: | || ::|  || :| || |:| |||| :|| |||:|||||::| | || | :|||||: || |:||:| || ||| ||||| ::|: ||| |||  :: ||:|| ::: : :|| || | |||::|:|: | |: | :|| ::|:::|||| :| ||:|:||  | ||    | |:| :| ||:   ||| || ||:  :  |::|::|| :|::|:  :||:| : |::|: |:|   |||:| |||:| ||   | :| ||  || |:|    |::|: | | :| : |:|: |||| ::| | ||||:   |:|| |:|||| |   || |:| | | | | |::|:| :|| |: |  || |||||| ||| || ||  |:| :| | |||  | ||| | : |||| ||||||: : :|  :||||| ::| ||| : |: |||  | | ||| | | :|: |||||  :|||  |  |||:  |||    : ||| ||||  |:|| ||| ::|:| | | ||:||| |||: ||: ::||: |||||:| ::: : |::| ||  |:|||: | ||:::: |||  |  |: | ||:  :| ||| |||:  | |::| |||:| |  |: ||: | ::| |||   |:  |:     || : |:|::|:| |||| :  | : :| ||||  |||:| || | | ::| |  |:| | :|:| | |||| :| |: |:||:|| ||||  ||| : | ||  ||||:| |:|||  |:||::||  :| | ||: ||:|  || ||| | |:|  |||  || ||: ||::||||| | | | |  :|: :|||| |:| |  ||:|:||: || || | :|: : | ||  || || | :|:: || |: | |   |||| | :||| |:| ||  ||| ||  || |::::|| || |||| :   |||  :|::|  |::||  |||| |::: |||:| |:|::: | |  ||:  | ||:  :  |: || :|||: ||:|  |||||||| |||: || |:|  | | |  ||: | :| ::::  ||||||  :|:|:  :|: |:| ||  |  ::|:||: ||| | |: || || || | | ||:||  || |||| || ||| || | || | ||  | ||:|:|: |  || | || |:| :|:||| ||  ||   |:: :| |||    || : ||:|||||:|  | ||| |:| ||  ||| : ||| ::|:|:: ||:||| ::||:|||| :  | :||: |||:||  |||| || :| || :| :::|||:| |:|::||  |: |  |||| ||| | ||:|::|   || |||  | ||||  :|: :||:| | ||:| |:|:|:|:|| :||:  |:| || ::|  || |:| :|| ||:|::| | |:|| |:| :|: | || |::|::: | ::|||| | |::|||  :|  | | ||| |||  |  |||:| ||: |  ||: :||  |:  :|||||||:|  | :| : |:|||| ||| ||:| | | ::| :| :  | |: :|:||:||  :|| |  |:| ||:||:|| |:| | | ||| |||::  | :|:|:| || |  | ||||| :||:|:  | :|||  | ||  |  :||| |  |:::| | |:|| |:|:||:| :|| |  | | ::|||  |||  | |   |||| | |: | :| |||  |:::|:||| || |:  :| || :||  | :||  ||  || |||: | || |:|| ||:||  | |   :|| | | | ||:|| |||:|  |:||    ||| | | :| ::||||:|: ||  : |:| ||:| |::|: |||  |||| || |:|| |: ||||||  | :| :  |:|||| |  ||:|::| |  :|:|  | |||:||: |::|| | | ||   |   ||| |||| |:|:||  :|: :::| : :::|| |||: | :| :|| |||::| |: || |:|:| ||:: ||| | ::||| || | :|:||  |:|   :|:|:| |||  :||  | ||| :| |:| | |  | ||||| ||| :|||:  ||::: |  | |:|| ::|: :| | | | |  |:|||||| | :| ||| ||| ||    || | |||:|:|||::||  |  |: | | :|:||| :| || |  :|::||: :||:||:   | |||:| |:  :|:| : : ::|||:|| : || :|| :::|:|| |:| || ||:: :| | | |: |:| || ||||:  |:| | :|| | :|| |: :|:| :||| |||: ||| ::|:|  ||| ||| ||  ||  | ||| |||| ||    || |:|:: :||   | ||| :||| ::||| |:  |  : | ||: || ||:| | :| :| : |:  : ||:||  || : ::| ||| |  :  |  |||| |:|  ||| || | ||  || :| |: |::|| ||| | ||  || ::::|| | :|   | |||   :|::| |  ||||| | || ||| || ::|| : |   |  | | :| ||||| ||:| | ||| :|:|: || : | |||||:| :|:| :|| ||| |  |:  |:|:|| ||::| ||:||: ||  |||:  ||:   ||::|  :|:|:|||: ||  |: |:|:||   || |  |||||:| |||  :||  | |  ||  :||||  :   | || | | |||||||  | | ||  ||:| :   ||| || || | ||||::| || | || |:||: :|::: |::: :||  |||| |:||:| |||: ::| :|| || :|: ||  |:|   ||| |  |||  ||||:|| :|:|: :|| |:|: |: |  : |:| || ||   :||: ||:| |: || | || |:|||  | |   ::||||||| |||: | ||  || | : | ||:::::|| : |: ||||  |||| ||| :|| ||:|| ||  ||| | :|::| :||  : ||:|::|| :| : || :|:| |     |: |||:||:: :|| ||:  ||||  ||  :| ::| || | ||:| :|||| | :| | ||: :|||||   ||:| | :||:||| :||| |||    | :  |:| | | || |:| || :|||  :|| ||| |  ||| || |:::| :| || : ||:| |||||| ||| : |||| | |: || :||:| |:|| |::|:  |:|||| | |||||| :|  |:  | ||| |::|  |  ||:|||||  :| | | :| || : | ||| :| | |||  |:| :  ||  | :|  ||||  :|||| || |  ||: |||| || ||: |:|:  ||:  |:| | ||| |  |||| ::||   | |||  |||| :|::|| | |:||| ::|:| || :   :|| | ||:|| | :|| | | | ||:||  ||:||:| |  | ||  || |  |::|::|:|  | |||  |  || |||  :|:: |||:|| || | :||   :| ::||   | || ||||  |: | | :|||  |  |||| | |||:||::|:||::|  ::||  | |::||   | ||:| ||||: |:::| ||||: ||| ||  :::| || |||  ||   |||| ||:|| :  |||  ::||: :| ::||||||  |: ||| |: : :|| | |:||| ||:|   | |||: ::||  : :| ||:|  || || || : |||:|| | ||:|| || | |||| |||| |: | ||:|  |  :|  |:| | |  ||| || |  |||  ||     | | | ::||||| |: |:|| |||: :|| :||:||||| : :   ||:|:|| || || :||:| |: |||| :|::::||  | ||| :| | |: | |||  | :  :|:||||  |||: ||  :| || ||| ||   | ||  |:||| || |:| ||| |||   |:|| | :|:: :| |||: : : : || ::| || ||    |:| | |:| :||||  |  | | | ||:| || |||| | :| || ||  |  |:||:::  | |:|:| || :||  | |:|||:  | ||: || |::|| :| || | |:|:||||||::| ||   ||: |  | |:: | ||||:   ||:|  |||:|  || |  || || ||:| | ||:: | |: |:::::| || | |   |::|||| |::::|| | |:|  || | ||| ||:|:| :||:|| | : |||  ||:|| ::| ||| |||  :|||:  :|  |:|  :|: ||    :|| ||| || :: |:|| | | | :|:||:   | :||| ||: |:|| ::| |: || | | ||  ||||:| |  |  |:|  ||| | ||:| ::|| || |::| | || ||:|: :  ||| |||:  || |:|||  || |:|  ||||:  |:|:||: :| || |:| |||:|  :|||| | |:|| || :|  |||:|| ||| |   :|: | |||: |:|:| :|:| : | || |:  :||| :|| |: |:|||| | |::| :| | |||| | : | :|| |:| | | |||| | |:||  ::| |::|  ||:|||| |||  ||:| ||     || :|| |  |  |:| | :   |:||||  :|   |:| |:||||||   | ||| ::||: ||| || | ||| |: | ||:|:| |||| ::| :| || |:|  ||:||  |:|||| ::|: |: |:| ||| |::| ||||:| || |  |||   :||::|| :  |||| : |  |: || ||||||||: | |:|| ::|:| | ||:|:  |||  | |||| :  ||| |:|:  |:|  |:|:|:|| ::||:   |||:| ||  ::| :|| || ||   || | :|||||: ||:| |: :| :|:|    ||||::||  :|| |||:|| |  : |||:|||  |  ||:| |: | || ||| | :| |: || | |  :||| :|::||| | ||: ::  | :| |: |||||| |:|:||| | || ||||:|| :|   : | ::|| | :: |: || ::||||  || |||: |: | |:||| :||  |:||| ::| ||:|| :::| |:|   | | ||| :|| | :| || |||| || || ||:| |::|||| :||: :|  || :| | :|||| :||    || |:|| ||:  |||: | ||  ::| ||:| || | | |: || || |  |:||||  ::|||||||  ||| ::| |||  || || :|| |  ||  ||||| ::|||  :|: || | |  :| | ||| :||  |   | ||  : |  ||: |||| ||::|:|:|  :|| |  |  |||  ||:|    ||||:||||  ||||:  | ||| ||:|  ||:|| :||| :|||| :|:: |: :||| |  |  | ||| |  |:: :|:| |  ||| :|| ||:||||  : | :|||| :|| | |  ||| || |||:::|    |:|||   :||| || ||| | | | |||||| ||  | ||||| |:|| :|| |:|: ||:||| |  || |:| :|||::||   |: :  |||| | | ||: |||:::| : | |:|:|:   |:|| ||  ||: | |  |||| |:||  :||     |:||:|: :|:|| |::|:  ||:|    | ||||||   | :|||:|| : |:  |||   || | | ||:| ||||:| | :|::|   ||: || |||||: || |||  | ||| ||  |||  |   |||| |||  |||  ||| : |:|::|::| ||| || : ||:||  |||  ||:|  :| | | ||  | ||  ||: |:|:|  |::|:| | :|:|||:  :|| || || | |||  :|||| | |:||  | |  |: |||| |||| ||: || |:| |  ::| ||||  | || | ||:|  :||:|:: :| ||:   | |  ||| | |  |:|:|: |||::||: ::|| | ||:|  | |  || |:: ||: :|  | ||| ||| |  ||||| |::||  ||  |||:: |:| || :| ||:::|:| ||| |||| ||:|::||      :|::|| |:| ||: :|||||:| | |:||:  |||  | | | |:::| ||  || | ||| |  ||||| | :| | ||| | | || || | : |||| || : :|| |||| : | ||| : |   | ::||| |:::| | || || |:|| :| || ||  || | :| | :||| ||||||: :|:::  :|| |: || :|||| |||||:| : : ::|| :|:| ||::|||:|:||   :||: | ||  | |:::||| ||:|||: || || :::| ::| :|:|  :|: || |  |||:  |::||: | || |||| ||||::|||| |  :||| ||:: |:|'
    'CCCCATCCGAAAGACTGACCACAGATAGTGCACTCCCCCTAAACTCGCCGTTC-A-AC-AACATCTATATTATCCAAAAATTTAAAGTTAAAACATACCATA-ATTCTATTCACAGAATAAATCCTTTCACGAGACG-TCCAC-ACGCGTTTTCTCATCCCAGCCCTACATACCCGGTAAAATACCAACC-TCTCAAATTCATAAGCGCACCCTCATCC-CTCCAAGGGCTGCCTGCACGGCAAGCGTCCACTCCATTC-CCTGGTCTCAAGAAGTCAGTCACT-AATAACTCAATGATCGCAATCGTTAAACG-CTGTCCGTAAACCCGCCTC-CTAAAATCTTCCCATACGTCTCATAGCGCACCTCAAC-CTGCAAT--ATCAGCACAAAGGAGCCGACGTCAATCAAAAAATCAACATGTTTCATAACTC-TCCCTTCAAACTCAT-TCCTGTACACCTTCTATCACAAGCGCCCTCCACTTATCTACTCATAATCACTAACCATCGTAATTCCACTGGGCATGGCGCATATATTCGGATCGCCGATACACTCCACA-AGTCCCATTAAACGCA-AAGA-CTTAG-GCTCAATTATAATAAGACCGTCTCCCCA-AGGG-AGATTCACAAATTTGCTTTACCTG-ATAGTTTCCCTCACCTAAATC-ACCCCCGA--GCTGAAGCCATCTGAGTTTCAAAGTTAATTGCTAGTATATAAAGCTA-CTACAGT-TCGAAAAAGT-CCTTTANTCAACTTACAAAACTA-AACCATAATATTAATGTACTGTTACCCAACATACATATCAGCGTTCCCCGCAGACC-TATCTTAATC-T-ACATCTCT-CATAACCACAGGTCATACTTCA-ATGAGATCACAGTATAG-TA-AAACGTTTAC-TTCAT-CTC-TAGTAAGAT-ATA-TAATGCATCCACGAATGCCCATTGTACTACCAGACAACATAAATTCCGGAAGGGGATGATTGTCC-CTTTCTACGCTTTCCTCCCGCTCTCAACCATCACAAATTAACTTAATTACCCCACAAGCGAAACTCGAGCCCGGA-ATACCGAAAGGCATC-AAGACCAAGCGCGTAAGGACCTCGGCCGAATCT-CGTTCAGT-CAACATACCCTAATACAAACTCCCAAATACATCTGTCATGCCTCTAGT-CACCGTCTGAATC-CAATAATTTTCCAGACTGGCTCAACCCCCCCCTCCCCGTCAATAGTTCATAGGCCACTCCCCATCTGTCGTAATAAACCATCAAAG-GAACAACATTCTCTT-AGACCCTATATGAA--ATC-ACATA-ACCCCTCCGTACGCCGCCATAACACTATAACCCGA-CTTTGGTAGAC----TCC-GCAA-TCAGTCTGT-ATATAATCAAATCCACACTCCCTCTGAACCAAACAGAACT-ATATTCACTAACA-TAAT-ATAAATAAATATTAACAAGCCCTACAGTGTAGTCCACTACTTTCTCA-GAAGGGTCTAACAGCTTACAGGACTAGACCCAATTATTCATAAACTCTCGCCCCGCAATTTCTTC-T-CGAGATAGTA--TAACG--AG-AAGAATAGAAACCACTATAAAAAGCTAATGAAACAACAAGTAAC-GTCCAAANACCTGTGACCCCACAACGCTGTTTTTTTCATCCCACCAAACTCTAAGTC--G--TCAGAAG-TC--TACCGTCAAAACAACCAGTCTAGCCCAACATT-ACTTATAAACC-CGACAATTTCACACTG-A-CA-CTCAGAAGCTT-ATACCGCCCCCACGTGTTCTAGCCTAAAGCT-ACCCCAC--CC--ACCCCTTAAC-CC-ATCAA--CGTACGACCAAACCCATATC-CTAAACATT-AGGAACGCTACCCGCAGAATGTCCCTATTTAGGAGCAAC--ATTCAAGACAA-ATTGCCACCCGACTCACTCGTAGAAAACCAAAAGCACCATACCAAATCAGAAACTGTGTTTC-CATGCATCTCAAACCCTATTTACAAGGGACTTCATAAGACCAAAG-ATCTTTCATCT-CGACCCC-AAAGCTT-GCAAAATCTCGCCAC-CGCT-CCCATATCCCAAAT-GG-ATTC-CTTGTACAGACATACGACCGTTATCAACAATTTTTATGCTT-TCTTAGCC---TCAT-TCACACTTG-ATAACTG-CATTCTCAT-ACCTACACG-C-AC-T--CCTACGGCGT-TACCAAGCTCGAGCGGAA-CTCGATACCAAATCCT-CAT-ACCA-GTGCTACTATGGCGATAAGCC-CACACCGACCGTCTTTGAAACTAGGTACAA-C-AATTTACTGTTCTC-TGA-TAGC-CCTTCCATAAG-AT-TTACCTACACAACTG-GTAAAAACGCAC-CACG-GTC-CAT-ACC-ACGCGGT-CCGTGCAAGGAAACC-CTACAT-ACAGACT-CCGAATAAATTAAACTCT-CGTAC-ATTTA--CA-CGG-T-TCGACTGTAATGTTATACTACGCTATTCAACACT-CATTAGAAACAC-TACCCAGCGAAGCTG-AGAGACCTCTTCTCATATGCTACCCGGAATCTTTACGTACCATACATATTTTATCAGCTCCTTGCAAATGCAACACAACGAAGAATCCGTCCCAGTTCATGACTA-CTAAATAG-TCC--TC-ACTTTGCATAAGTTCAACTTATG-ATCTC-CAAGTATTC--TTC-ACCCGAT-TGAATATCAAC-CATC-TCAAGTATAAACCGTTTA-AC-C-CTTA-ACCCCCACCACGGCCGACACAGATACTACAGA-AGGCACAT-AAACACC-TTATC-GTTCATTTCACTGCCAC-AGAAACAAATATCGAATACAATTATGGAC-CACCGTAAAACAAA-ACAACATAAGCCAAATGTTCC--ACCATCATGAAAAAATACTAACAGCC-TCA-CTACACAGGCCAGCAGGACTCCCCAATTCCTCATATA-ACTTAAAACTGCCATGCGTGACAAAATACCTCATC-C-ACAAACCCACCAAATATAGATTCCGTATTTCCATTCAG-TCAGAGTTACTCAAACCCAGCAAACACCAGCCCTAATTAAAAAAAATAGC-CGACGCGCACATCCAACCCG-GATGC--TCCAGTGATCACGTAGTTTCTACACAAGTGGCGAACTTACGTAACCTACCCATACTAAAAAATAAGAAGCA-C-A-ACCTCCTGTATCCGCAACAGTATTGTTTGCACTGGGCACTNAAT-ACTTCG-AAGATAACATA-C-AATACTCGGACATTCTAGC-T-CC-CCATCCGCGATTACTAATCGCTCTTCTATATTCTAACTGGAACCCGTCTGTACCCGCCTATTCACGGACTACCCCATACGT-CGGGTATAGCATCAGCTAT-TC-CTTCCAATAATTCCACCGCGC-TCT-CCTTAGGCCACCGCTGGTCAACAC-T-C-TTGCC-CCACAACTTACAAG-GA-CCATGGCCTAAAAG-CCCGAAATGCGC-GACTCCTAGCGTAGGAGGTAC--GAGTCGAAACTCAATCATCTAGCGCCC---A-GCTTTCTCATG--AGCTCAAATGTTA-C--G-ATTC-A--CCTCACATTACAGATA-G-AACACACTGC-CA-C--CTTAGAAGCTTCC-AACGCTCAGACAGCTATTCAGCTTTCAACTA-TGCT-GCTGCACCATGTCTACC-GGAAACCTTA-T-A-CT-C-AAGCACGA-TTCCGTAATCCG---C-G-TACAATCACCCCGAC-CATCCTC-C-CAGTTACTGCCCTGTCTACGGAACTGTTCTTTGCTAATAG-C-TCCCAATGGGTATTATGTCAA-TCTTCAGTAAACTTACTC-A-CACT-TGT---TTAC-C--TCCCCAAATGCAACA-ATCC-TTCACA--ACA-CGC--ATCTCTCGA-CCAATCGTAGACCAAGTACCTA-CTCCC---AGTG-TTA-CCCTTGAACCTTCCATTGTGGTCCT-GACCA-AACACG-CCAC-CCCA-CAAC-CCAACGCATGCCATAATTGAGACGCC--CAATGAC-AATATCATACCATC-CACTCAAATTA-CAGAAGTAG-TCACA-A-TCTAG-CAA--GCTGGGTTGGCAAGTGCCACCCTG-CGCAAAAACACATGTAAGTTTGTATCCTATCCCACATTTCATCTGACCTCACACCAAATGTACTCAACATCTTTTAAGACATTCTGCCCCCCGTACAACCCTATCAAGCGGTATATACAGAATGCAAGCC-TTCAATCTTATCTTCCGAACGCACAACAACCTTCCGTATCTACC-ACCAGTCACCAGTGTATAGAACCCGATTATAGGTAC-GCCTCACTACATTATGCTCCCCCATCAGCCGCCCCCAACGATTCTCACCCCACGTTTTAAACGTCTT-AACAAAGCTTCCATCACCACCGTACTTATCTCCACCTAATACGGCCCCTCAGACACGACGCCAGACACTCTCACGGGCAATTGTCTCAGCATTATGC-TCTCTTATACCGTCTCTCACATCGAATTAATCTAGTTCACC-ACTCTCTTCGG-ATCTTGAAAATCGACTAAG-TAAGAACACA-GACT-TA--ACCGA-AATCGCCTTGGTGCAGCAGCTTTCCTNCGCAGAAATTGCTC-TTCCGTAGCACTAACGAAGTCCTTTTCTCCATCTCCATG-CATCTCTGACAAACAAAACACCATGACCC-CCCAAATACCCCCCTAAAACTGATTCCCACCCCA-C-ACAACATACTGAATCGGGAAACTCACACCACGCCTAATATCCGCCTAAGCATTCTGTTGTA-CTCAGTAGGACTGTGTTAGCAACAATAGAAGCGACAATATACGACCCCGCTCGAATTGC-CA-TA-CCTTAC-CTAACCACCTAGATCACCAATCCA-GCC--TACCACCAAGCTCGTAC--TCCCGAATACTAGTTTCCTATAGGGTGAAATATGTC-CCACCATCAATACAGCTCGAATTTGTGCCTGTCGTCACACTTGCTATGCTTTA-TGCTTAACCATCAGTTCG-CAGACTTTCACGACATAACTTCACTTACAGAGAAAAACGTCCCTGAACAGT-C-TCGATAAGCACAAGGACACCCGCCGCCCT-C-AAAACCCCCGTC-GCATTATAGTCCCCTGCTTCTTATACACATGGCCCCCCTCCACCTTCATAGTCTAACNCCGCCCAAACATCTTCTGTCTTACCCTGGTACGACAATCCTAGTTACTGAAAACTCCACCC-CTACATAAATTGCATC-CCCGCGCCCACAAACACGCTCCTATACAAAGCA-AAATACACGATGGGCCTATTAGC-TCTCTACTCTCTTCAACAAAGACCGATTACCCCAATCAAGC-CTCCACTCCTGGAGGCCCAAAG---ACT-CTTGAAATCTGATAGTATAAATTACCA-ATTCTACAACCTTTTTCAATTTGACCCTCA-TCCATCTAAAGCTATCCACGTTCAGTGACACCCGCCCAGCTTCCTGATAC-CCTTACCGCATCTTGTTCACGTAGGTTCCTTATAGGTTCCGGACTGCAA-CGCTTTTTTCTCATCCATGCACGCCACC-CATCGGCCACTGACAATTGTGATACGCTCCATATAC-ATCAATACTAAACGACACCTCCGTAAACCTA-CCC--CTAACAAAAT-TC--CATC-CACA--CAATAATCGTC-TCGCTATTTC-A-A-CCCCTTCACCC-CGCACGTTAGACCCTCCCCTGACCTCCTATATTAAGGAAGCACTACTACA-ATGCCTAGACTCCCCCCTGTCACCTTCT-ACTCGCCAATCACCACT--CGTAGCATACAGAA-TCCTCGTTTAACTATTGCTACACCAAGTCCGACA-TCTACACCGATGGCTTTCTCCCCCATTTTAGCAATCCGTACCGCCGATAACA-AACTACGTCACCAGTAGATCAGAC-ACGGGAAACCATCCCATCTATAAGTAACCATCCCACCCCGATTAA-ATAATGCGTTCATACAATGCATCCCAGCGACGACCTGA-ACGTCTACGTTACTACGCTTAAATGCACTTCAGTGCCCGGCA-ACCACTAACACGTGCATCG-CACTATTGTTTTCTGACGCAATTACAGAACGCACTTAGACGTTTTATATCACTACTACC-TATTGCTTACTGCATTAACCTCTTCTTCCGTTCGGTAAATCACACGTTCCATCAA-GCAAAGAAATCTAACTAACACT-C-CATGCC-GGTACCCTTCTTTGCCATGTCCTTACTCCCCAATG-CACCATCTCACCCTTCCCTCCTTAAAAACCGTTGAA-TTTAAGCTTCCCCCG-AACTCCCTAGTCATCCCCGCCA-GCAGCACTGAAAATTACTATTAAACACATAACCGCG--AAAC-GAT-A-CAATC-CGCCAATGTAGCAAACTAGGAACACCTACCCGTTGATGTCTAACCAACCAGCTAGTTGCGATTTCTTCCCCCGAAGAACCTAGAAAGCGAACATTT-ATTTTACAGAAAGTACTTTC-AGCAGCCCAAGAGTACCCTTCTCGTACGAAATCTACAAATAAAATTTCT-TGTACCTAAGAT-TTA-CCCTT-CTGTAGAATAACACCACCCCATAAACCTTTTCCAGAATCGGTATATTTTTAATATTTGCCTATGCCCAAGTCTCCACAAGAAT-CTTCACGAAGCCTCTACCTATAAGTAGTTATATTTACCGACCTATTCCCCCAATACCATCCTATACTACTTCAATTCAAAGCAGTACTTGACCAGTACCC-GACAGTCTTGCCGTGAATCATCGATG-CTCT-ATACCCTCACAGTCAGAACCATACTG-GCCATCAGACCC-GAT-AAATCGCATCATCCGCGTAACCCACCAATGGCAATGTACATA-TC-T-CCGTTCGATCCT-C-ACCTAA-C-GAA-A-TTATTTCCATTGAATAGA-GAGGATATATGTCCATACAACCCCTATGAAAGATTACAGACGCATGTTGTACCAACTTACCTGCTTGCACAACCACACTCAAACCCAGTCTCATGTTCCACAAGTTTACCTCAGACCTTGCTATAGGATTCATTGCC-AAACAACGCACATACTATAT-TCTTTT-TCAAACTCTTATCCATTCACTTCATATTACTCCAATGNC-TCCTTTG-ACAGACAGGTCCACTTACATTCTTATCAACCCT-GTTGT-AGCTGACGAACTAGGCTTACTC-CGAAATCTACCCAATCATATAACCCTTC--CCCTGGTGTTACATACTTAAC-CTACAGAAACAAATCTTCAACGCGTTGGGCTACCCGCCAATTTCTGCGAACGTACTGCATAAAACCTTACCGCGGAA--CAAT-TACCGACACAATTTAACAGCCGGCCTC-CGTTAAT-ACAACCTTATCAATAGTAA--GCGTTACCAAGATTCTCTTT-CCT---CTA-A-CCCAAACACA-CTT-TGC-TTAGT-AATCTTTTCGCCCACCAATTATAC-CTTATTCCCAAACAA-TGGGA-ACAATAAGC-ACCTAGCTCAGTTAACT-AAAC--CAGCTTAACTAA-CCTCT-TCAGTACTGCTAAACC-CTACGACAAACCATAC-T-TGCACAACCCACTGACC-GTTCTCATC-GATCAGTTTC-TGT-CCTCCGAATACCANCAT-TTCTGGACATTACCTACATGTCATAGGGAACG-GGAACAAGCGCACTTACCCC-CAT-C-TAATC-CCG-CCCATTACAACAGAAGAATTAACTCCCACGCTTGAAGTCATCCGTGCACTCATAGTC-C--CAC-G-CAACC-CATG-GAAACCCCGAACCGCTA-CAACTTGC-CCTGTTCGTACAGTTCCATTACTA-CTATCC-TCTGACTATAAT-ACAATCGTTTCTCACCTTTGATC-CTTCTGTATCTTGCAGACAC--CCCCCTCAAAACGA-AACAAATCTGGAGCATATACCTGAAACACCGGACATCTTAAAGCCTCGCA-AAAC---CAAAACACTCAATT-TATTCATCCA-CATTC-CTGATAG-AAACATTAGAATGAGACTCAACATCGCCAGC-TTC-GTGGCTACG-AACCCACCCCGCACTCAAAAACTAACTTATTGAATTCAATCTACATCTACTTTTGAGGCTACTCGATCGGCCAGATAAATTATTACGC-ACCATTCGCACTAAAC-CAGCCTACACCCCAATTCGCTATTTAAAAAGTTCAATTGCACC--TCCGTGATTTTCACCTTCCCCCTATCAGTCAATTCC-ACATGAGTCTCAAATCCCCGATCTACACGCAACCTATTACCCCCTACTACTGCATTTGTTAAGCCGTTGTATCCAACAGTAAATTTACCGTAACTACTTCCTCTGCCCAAGCAG-CCCTACAAAAACTTATACACCACCCCTTTAGTCGTGTAACTCGGTCTCCACCCGTCA-AAACTACATACAAATGACGAGTAGGTAG-CTATATCATAGTATAACATTCCTTACATATACGTTTCAACAAGTCCACATCCTTATCACAAAGTACTATAAATTCCAGAAACTAAAGGAACTCACGATCCCGTAGTCTGAGTAACAAACCGACTTTTCACGACCG-TA-CTAGCCATTCTATCTGTACATGCTCCATAT-GT-GTCTTTGC-CTG-GAACAGCCAATTGCTCCCCCCGATTCGTCATAGC-GCAACACAACA-ATACTAC-GAT-GTACCCCAACTGCTCACGGTCCTATCA-ACACTATTTAGTGC-ATCTACCGCTTCTACAAAGACTAACATCA-ACAACATACCCAA-TACTTACT-CTCAAAAT-AGTC--CGCTTCCATA-CTACTCAACCAAAACATACTGACTTTA-ACCTCCCCAA-CGCTAAAAAGCCCCCAATAGACCACATTTTCTCTCCCAA-ATCGTCCNAA-CCGAATTGAACACATTTCATCAAACCTCTA-CCGGGATCACAAATACGGTATTAACCGC-CGGTCCCTCAACCCTCTACCCTA-ATAGAACAATAATCGCCCTAAGGAC-CTCATAGTTCAACCACATTGGTCCACACTACTCACATTAATTTCG-CACCTTCTACAAGTCTTGAATC-TGAGCTCGACCT-CAGACTT-TATGCAGTGTCC-T-CGACTG-CGTTTATCACGTCACCTCTGCACCCAGGGGAACC-GCT-ACACCTATTCGGACCCCCAAACCTCAAGC-TCCAACTCTCACACCGTCACTGT-GCG-ATAGCGCTA-AGCACACGACAATACAGTACTGTCTCCTATATGTGATCATTATACTAGCCAACAAAAAATCCACCCTACATTTCGACTGATCA-TATCATTATACCCTATACTCCGTCAACTAAAGTAACACTGCTGTTCGACAACACGAAT-TAT-ACAGTCCGTTTCTA-TCACCAATATCCCAGACCACTGCGCCAAAACACAGGATACTGACCTTGGCAA---AAGGAGTTT-CCCAGCTACG-AACTAACACAAGTA-AGCAGCG--ATGGCCTCTTCACACTCCAACAACGCAACTACTCCATAGCCCC-CAT-TCCGT-ACATCTAAAGATTACCAGGACCCCACTAAACCAT-A-GTCCCATTCACAG-GCA--AGCCACCTG-CTT-CCG-T-T-C-CCC-CAAACTC-G-CCCCTTAGTACC-TGCGCCTATTATTCCTGT-TGAGTCGTC---TCCTTATAAGCTGCTGTCTAGGTTAAAAGGCAC-A--G-AAAACAAATTCA-CGAGACACACTGCTACAATTACCTGATAATACAAACTTAC-CTTTCATGCTC-CCCCTCGGTCTTCAGAACTAAATCTAGTAC-ACCCTGTCGCAAATTCGCTTAGCCAGCCGCTTAGCACG--GCTTAC---GGCTC-CCACATTTATATATATTTTAAAAATCCC-T--CCAAAAATATAAAAG-GA-CGTTTATTGAGAATAGCCCCGTCACC-CTTAGTG--GAAAC-AGTACTTAAC-CAACTCACGTGCTAGTACCTACTTATGTATTGTC-ACCTTAACTAATAACCTAG-CATACAACGTCACTA-CGACC--CTC-CTAGATTCAC--CCACGAGCACT-CG-C-AACGCTCATAGCCCAAAGC-CAAACAAC-AAAGCCTCT-C-CGGTGTAAAGTGCAT-CGAATCACATAAACAATG-ACCCAGCCTATA-CACTACACCGCCCGTCCATAACC-T-CTCG-CCCACT-GCCTTACTACTAAGACAAC-CACACCAGCACCGGGGTACGAGAGCGATACCTGACGTCCAAACAATATAAACCATAATAGGTTCGCTGGCCCTATGCCATACT-GCACTACTACAATGTCTCTCAAGAAATTAATCAGTTCCAAGCAGATTTGCNCCAGTTACTACCACTACCTGCCACGAATAGTAGAGAGGATACAAACATTTCCGCAAGACTGAAATGGCCAAAGCCTTCAGTC--CCCATACGATATCTCTGCTGTACACACGAAATGTGGAACAATCATTTACCTCACAGCGAAATTCTCATACCAGCAAAACAACCAGCCTTTCCAAGAACCTACAGCTACCGGCTCACTAAGAACATTGTGGAGCTTCAAGAG-AC-AGCAGAGTCAATTGTAAGTGAATAAGCTCTACGCAATGTGTAAAAATTCAACGTCTATCCCAATTGT-CGTA-CATACATATTTTTCTATACATTTCAAGAAATCACCCCACGTCCGTCGCGCCCCGGATCGTAACCGCCCATAAATCAAGTTCAAGGTCTACACAAGCTCCACGAGATAATCATGGAGGTAAATTTCGCACGACTCAATTAAAAACATACCTTTTGCCTGT-CGTT-TTTCT-TCAAGTATAACCAATTTAATG--TAA-ATCAC-C-C-A-ATTTGGCACCACACGCCTC-CCCCACGATAGACAATGCCGCNTCATGATATTATCCT-GCACTCATCACATATTTTTTGATACTCAA-AATCCTAAACATCGTCCCACCCTTTATTATTGTCACAAATGTCCACCCCCGTATGCTCCTTCCACTTCATAAACCTACACGGCGGTCTGCTTTAACGCAATCACACAATGAGCCGATAATAACCATTCTCATTCACCC-ATTTCT-G-CTGGATCAACCCATCTTTCGCAATTAC-CCC-CCGATTTGTGTGTGCATTGATTCTCCCAGCCATCCTATCTCAACCATCCCCCCTCTCGTATGCAAAAACCTT--CTTCACGATTG--CACACTAC-AGATTTA-TCCTAACCAAGGCGC-CGCGATCACACAACAATCGCTCA--CCCGAGTTTC-CAAACTC-TACCC-CATCCGACC--GC-ATAAGCTTACTCAC-TTT-A-CGC--AAT-CTC-AAGCACTTCCATGTACTGC-ACTAAACTTCCCTCACACGAATAATCCACTACCG-TAACAGG-TCATTACAAGAAAGTCATAAC-CAGAGCTGTTAAATTATACT-CTGACCCTTAC-CTAAAACTATTACC-TATAAGTCACCCTACCCGGATCAAGTTGGTTTACATAAATACCAACGCATAAACACTGAACAAGATACAAAAGTATCTAAGTATAGTGAATTACAAGCAT-AGATCACAAATTACGACACCCTACCCTCCAAATCCGAACCCCCGCGGACTGAAGAGATACAGTTTTTCATTACCCTGTGCACCCTAC-GC-C-T-GAACAG-CCTGATATCA-ACA-CATTATGTC-CTTCT-CCTATGTCTTTCCCACACCAACTAGACCCCAACCACCGGAAAGGGAATCGGTGAATAATCACTAACACTACACAACAAT-GCTGTTACATTAATATCACTTCTATCCACATTAACAATCGCCACTCAACGCCTTATCAAACGTACAAGACAAC-AAACAAGAACAAGTACAATACATTATCC-CGCTCGACGAAGACGCCCCGTCAGAGTCCATAAAGTTAGACCACCAAACTTAAAACAACGGCAACCCTGCCAAAGC-CCACAGATTACGCAATTACAC-TACA-CACGTGTCACTTACTACAGACTGTTAAAAGGCAAACCC-CACTTAATCACTCAAA-AAGGTGTAT-ATACTCCACATACT--AAG-GCTTAACATCC-A-TCTTTATCTCTACCA-AAAG-CATAAACGAGCTCCAAAACACTAGACCTGAAAACGGAACACATCGGAG-TAC-TAACCACGTCTAATATTTGTAAC-TC-A-CCCGCATTCAAAACATAC-ACATC-CAACCCGATATCT-TAAAAAAGCCCCTATCCGCTCATC-TCATCTTTCGTTTCTTTAAATGTC-T-CCCCCT-G-GATCAATATACTACAAATNATAAGATATCCATCAA-ACAGGGAACACTAGATCCTCCCAGCTATTCACCCACTTCCCTGATATTAGCCGATATCCTCATCTCGCTTGATCTTAGGCGTAGTAACGTTCCCAACAGATAATTCAACACACCATGAT-TTTCCCAAACTTAAAATCAATACTTCCATGGAACACCCTAGAACCTAGTAGATACTGGACTTCCA-TTAGTAAACTAGCACCACTTTCCACCCACTAACTCTAAAGCCGTGATTTTTCTAGAAATACCCCATGGAATTCCGGTTGTACCACCCCAAGACCACAGCAATAC-GTCATCAATCCAAAACGCCACCCCTAATAC-TTCATGCATCGCCACTACATAATGCATTACGATTGATAACATACATCAGGAACCTCTGGACCTAAGACTTCGCCTTACTTAGCAGAGCGCAAATGTC-CCCGTGTTTTTCTCTAAAACTCCCTTAAC-ATCCGCACGTGTTATAGTACTTTCCGACTAATACATCGATTATTCTACCAACGTCGACAA-AGGCGGCCTTTCCGCACGGTTCTTTCATGTTATCTATTGTTCAACGCCCCTTCAATCCTTAATAACATGTCGTAACT-ACTGATTTCAT-ACAA-CAAGTGGCCAAAGAGTAATACCCTCTGCCTAGTGAGATATGAGCCAATACACTTGAACGGGTATATCAGTTTCCAAAGCCACAGCCGTTTACTCGGACTCTTCTCCTCGTGCTCCGCCTTGCTGACCTGCGCTACCAACCCCCTCCTCTTCACACTATTATCATTAT-ACTGCA-AATAACAAGCTAGTGACCGA--CGTTTTACCGCAACTAAC-GTCC-CAATTCTCTCAACTTGATAATCTC-TCTGATATAAATTATGTTTACGTCCGACCGACAACA-CCC-AANCTGGTAATCGAACTTTCTACTA-ACCTGACACTAAATTGACTTA-CATGTCCAGTATCAGCCCTTATACCTTCTCATGCTAATAATCCTTCCAGTCAGATCGCC-C-CT-ACCTCATTGATGACCTTACCTACG-ATCTGGTAGCAAATTTTAAATCCAGCCCCCGGTCACATTATTGCAACACCCCGACCGATGAGAGAATGTA-ACGT-TAACACATCTAGAAACGTTCNCCCCACCACTATCGGTTACAGGAA-AGATAATTATGAGCCCACATACAATAANCCCGGAAACAAGAATGCACTCTCCTCAAAAGACCGCCTG-AACTCC--ACATTAGTTGCATCAGA-CC-AGTCC--A-G-CCTTTGCATACGCGAT-C-TG-CCTGCCACGAACTAATTCAGATAACATAAATCTAT-TAA-AAAATTACC-TAC-TAT-T-AGAACCACTCCGCCCTTTGCTCAC-TCGC-AGAT-C-CCCCCTTTTTTAAATGACATACACCTCCATCTT--TACGTC-CGTACAGTAAAATTCCCCTCCCAAAATGCGACCAGGCCACACCAGAACATAATTCACACATTTATGCGCTCTGACTACACTTCTTCACAACCAAGGACCTGTCTTAGTCCAAATTCATGTCCTCACAA-CTTCTT-ATATAGAAATCA-AACTCCTCCAACTTAACATTCGCCAAACAATTACAAACAGAC-CCGTGC-CTTACGTCTGATAAATATTAC--TCACAGCNCAAACGCTATCA-AACACATTCACGCAAAGCC--G-AC-CCTATG-AAAGTGTTAACCGAACGCCCGAC-CCAGATCAGACATTCCAACGTAAAACTTAAACCTTCCTATCGCCT--TCCCCATAGCC-C-GTATA-GTTAC-TAGAGATTAGTCTCGTCACA-TCATGACTTAATGGCATCCTGCACCCCCACATTCCA-C--CATGAATTACCT-ACAC-GTCTCATTAGTGTAATAACATG-CCTCGTCA-TAACATTTTTGATAGCAAATTTTCA-ATTTGTATG-GCCCACGCATCCATCAAAGATAATTATTATTTCTCATGACTGTATCTTCATAAATCTGTTCACCCTGCCACGGCCTTATC-TCGTATATCCGAAACCA-CTGCA-GCTTGGCCATACTCACTTTATGCAACCCA-TCCTGGGACACGCCCTACGCGTGGCATGCTCCATCAAACGTGCGGCATACCTAC-CTCACCTCCACACCTC-CCTCGCCCGGGACAT-ACTT-GG-GTAAAACCCCCGCGTTTGCCCCCCACA-AC-A-CACC--TG--CCCA-A-A-AACTCAATTACTACT-GAAACTCAGGATGTCGCCCACTTGATAACACAGTATCTTAACACAGTCACCAGAGATAAATCCATAAATATTCCCTCCTCAACCT-TCGAATACCG--GTAAGCTTCTATTCNCCC-GCAGTAGCAAGTTATATACCAACGGCACATTGCGATCCCCCAAAAAACGCCTACTCATATGTGTGT-CC-GC----GTC-CCCAATAAAACC-A-CATG-CTAAAAACC-TTCT-G-TCCGA-ATATAAGTATCTNGGTGTCTTTTATC-GCGAGCTAACACGCCCTGAAGAG-TAACAGAC-ACCA-GAAAACTATCAACACAGGCCCAACCCG-TAGCCTACCTG--CTCAAC-CCGAAAATACA--TTTCA-AATC-T-AGGTCTGCTCCGTAAGCTTCCG-A---CGGGCGC-GTTTGCACG-TTGCTAACCCTACGCCTGAAATTCGCTATCGCACACCGCGTAAAAC-AACACTACGGGCCT-CGCCGCAG-CGCTTACGGAGATCTAACGGTCAGGCTAC-GCCTCTATTACACCCCCACTTCTCACATACGATCACTTCCGC-CTTAGGCACTTACACCCTCTAATTTATCACATGAAATACAGCCATAAGGCCACACTCTTTTCCCTAGCGATCTTCAGCTT-CCAGCGC-CCACCACTCTAACTGCGTCTCGATGTTCCTATTACAGTTA-TATCTTCTACCTTCAATA-TAA-CTGTTCTGAGAACAGAACTCAATGC-TCAGAACCACGAGTC-GAT-CAGTATACT-A-TAATAAATC-TAAATTAAAT-TCCCATGCTGA-CCCA-CTTAGGACGTCCCAATCCATGCAATCAACTTCCTC-AA-CC-CGACAAC-TTTC--CAT-TCTGAGCGCCGCCCGTCCGCTACATAAATCC---CATAGTTTACACATACTCCT-AAC-TA-CTGATCACATCTTTATG-TCACC-AAT--CAGAA--AGACACGTTCGCAT-TAC-TAT-AAAG-ATTCTG-CA-CTGTTC-CCTCCGTGT-TCCAACCTG-TCACG-CATGGAAACACTCTCATACCAAGGCTTAACTAGTATTC--A-TTGTAAAGATA-ATGGCCCAATGCTTTCTCC-TTCAAATTTCGTCTACGTTGCTGATTCC-AC---CC-TCAAACCCGACTTACCATGCTTCTA-T-TACACAACG--AAC-CTTAAATCGGACCAAACGC-CT-CATTCAGT-CCGCGACATCTAAGTTCCATATACTTAACTTGTTGAGACCACAAAATAGAAC-CAGGTCCCCTTTAAATTAACAACCCCTATCCACCGCT-CGATGCCTAATAGTGCCTTAAATG--TTCTTAG-CACCCTATCAGAAC--TAAACGATCTAATTAAAGTATCATC-CCAGAACCAACGAACATTAC-CCCTCACCCATACACCAGAACAT-TCCAGATCCAGAAAAC-ACTCATGAGCGTTAGATGCC-CCACCATAA-TCGCT-AAC-CCACACG-AAATAAATC-CTGT--CACATATCATATCAAAGCCACCTAAGCCCTAGCGTCC-TACCCGTACCCCGGAC-AAAAGCCAAAATCACACTGGCCCTGCTTTCCGAAATTGTCAGTATG-AAAATTTCACCAGTATCACGTTTTGCAGATGCCCGAATATCCCCCAAATTAGTCCAACCGCCTCCCTGAAAC-CAAGCA-ACTACAGTGCATATATCG-GTTTAATCACCCAGACGCTACTACGTAATAGCCAA-GACCTCTCCTTGCTTCGACCCACTAACCCA--CTCTTC-CTTTTAG---CGAATC---GACA-ATTTCCTA-G-CACCAACCCCCTCACAATCT-CACCTGAACAACTACCAACCTTGGACGTAATTAAAAACTTTCGTATCCATACAGGCCAATTGAAAAACAATGCCTCGCACCCGTCACTAGCCAACACTCCATCC-CGCAGGGCCGTAGGCATCACTCACATACTTTTTCCGTGCATTTGCTTTACAGGGCAATCCGCATTCAATTAGCGAGCCCATCAATAT-TAATAGAATATTACGGCTCTGCTATTCCTTCATAGCTCCCCCAGCTAGCTACTCTCCCCGACCATTAAACAGACCAGGCTGAAACTTTTGAAACTAAACCACAGTACGAGCTACAGAAGTGACTCTGAATCGGCCAGTGGACTACTTGTTAGTGCGCCTTCCGGA-CAACCTACTACCCACACATTCCATAGGAGAGGCTTATACTCACCTATATATGGAATC-CAACTCCCACCTCCCTAAACAC-CAACGAAATAT-ACTCCCACACCATGCACG-GAATGATGTTAGTCTGGAATC-C-CTTCC-T-C--TA-CTATCTACAAAAC-AGT-GCA-A--ACTAA-CAAGTC-CCCCCCCAA-CGC-CTAAATTCAAAAAC-AAC-AACAGCCGTTACAAGGGTAATTAACCCAGTTAAAATCCGCACTCGAGACAGCCCATTATCATGCTAGACTTCCGCACCAGACAGGACAATCTCGCATACTCCTGTGGG-AACTTATCTCTATAGGGGTC-CCATTAGC-TCTCTCACCCCCAACATTAATATCTGAGTTAGCATAAACACTGGTCCTACCCATAAAACCG-GACCCCGCNGAGACCAACGTTCCCAAACAATAAACTTCTATAAAGTCCATTT-GCATGGAGTTCTTCCAACGCGTTTCCATTTGTAAC-CAAACTACAAGTTAAAAC-AGCAATGCAAATATGCAT--GACTCCCTTCAATGCACTTGAACGTTGTACCACTACCAATTCCCAACCAC'

>> showalignment(rlocalAlignment)
>> 
>> 
>> [rscore, rlocalAlignment] = swalign(hippoPORF_random, humanPORF_random)

rscore =

    23


rlocalAlignment =

  3×84 char array

    'LHITITSLTISGLLIGLTPVTPV-EHILLSTFPLHANNAILPMYTPIISAILLLLMPFLLVSIIQLTIAIQLINRSLNIKLSII'
    '|:: | : | |:||     |:|:  :  |::   : |:: :|:|| |::::|:| :     : ::  :|| |:: : |  | |:'
    'LNFEIET-TYSALLQLYFLVSPILAYPTLTS-GANNNGTKFPIYT-ILGSLLFLALSQKTSTDLNTLLAITLLTLQ-NTPLIIV'

>> showalignment(rlocalAlignment)



codoncount(hippo)
AAA - 197     AAC - 175     AAG -  75     AAT - 153     
ACA - 159     ACC - 140     ACG -  47     ACT - 131     
AGA -  63     AGC - 121     AGG -  64     AGT -  63     
ATA - 157     ATC - 151     ATG -  67     ATT - 106     
CAA - 170     CAC - 142     CAG -  74     CAT - 153     
CCA - 116     CCC - 140     CCG -  42     CCT - 142     
CGA -  32     CGC -  37     CGG -  30     CGT -  39     
CTA - 192     CTC - 121     CTG -  52     CTT -  96     
GAA -  73     GAC -  53     GAG -  42     GAT -  48     
GCA -  80     GCC -  80     GCG -  21     GCT -  52     
GGA -  44     GGC -  47     GGG -  28     GGT -  23     
GTA -  52     GTC -  36     GTG -  29     GTT -  42     
TAA - 113     TAC - 111     TAG -  82     TAT - 106     
TCA - 114     TCC - 114     TCG -  38     TCT -  94     
TGA -  72     TGC -  44     TGG -  30     TGT -  46     
TTA -  93     TTC - 100     TTG -  41     TTT -  72     

codoncount(human)
AAA - 167     AAC - 168     AAG -  72     AAT - 133     
ACA - 134     ACC - 184     ACG -  42     ACT - 157     
AGA -  54     AGC -  95     AGG -  60     AGT -  58     
ATA - 119     ATC - 132     ATG -  52     ATT - 109     
CAA - 137     CAC - 154     CAG -  68     CAT - 150     
CCA - 145     CCC - 198     CCG -  52     CCT - 186     
CGA -  43     CGC -  55     CGG -  26     CGT -  26     
CTA - 165     CTC - 140     CTG -  74     CTT -  93     
GAA -  61     GAC -  57     GAG -  49     GAT -  44     
GCA -  81     GCC -  95     GCG -  15     GCT -  58     
GGA -  35     GGC -  42     GGG -  21     GGT -  24     
GTA -  39     GTC -  26     GTG -  13     GTT -  38     
TAA - 165     TAC - 118     TAG -  96     TAT - 104     
TCA - 126     TCC - 114     TCG -  42     TCT - 101     
TGA -  64     TGC -  44     TGG -  35     TGT -  32     
TTA - 100     TTC - 108     TTG -  49     TTT -  74   

orfs_valid = seqshoworfs(hippo, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 94);

orfs_valid.Start

ans =

        7945       10186       12064


ans =

        2750        3965        5339


ans =

        7032        8625       14169


ans =

        2299


ans =

     []


ans =

     []

orfs_valid.Stop

ans =

        8623       11578       13582


ans =

        3704        4955        6899


ans =

        7713        9408       15306


ans =

        2833


ans =

     []


ans =

     []
     
orfs_human = seqshoworfs(human, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3], 'minimumlength', 78);
orfs_human.Start

ans =

        3316        8530       12346       16570


ans =

        7589       10763       14750


ans =

        4581        5907        9210       10326


ans =

        1897


ans =

     []


ans =

     []

orfs_human.Stop

ans =

        4264        9208       14149


ans =

        8270       12185       15908


ans =

        5514        7446       10056       10767


ans =

        2419


ans =

     []


ans =

     []

seqdotplot(humanProtein, hippoProtein, 3, 2)

Warning: Match matrix has more points than available screen pixels.
         Scaling image by factors of 5 in X and 9 in Y 
> In seqdotplot (line 261) 
hippo_orf1 = hippo(orfs_valid(1).Start(1) : orfs_valid(1).Stop(1))

hippo_orf1 =

ATGAACGAAAATCTATTCGCCTCTTTCATTACCCCTACAATCCTAGGCCTACCCCTAGTCACCCTAATCATTATGTTCCCAAGCATACTATTCCCAGCACCCACCCGTCTAATTACTAATCGCTTAGTCTCCATTCAACAATGACTAATCCAGCTCGTATCAAAACAAATAATGAACATCCACAACCACAAAGGACAAACTTGAACACTAATATTAATATCCCTTATCCTATTCATCGGCTCAACAAACCTCCTGGGACTCCTGCCACACTCATTCACACCCACCACACAACTCTCAATAAACTTAGGCATAGCCATCCCCCTGTGAGCAGGCACTGTAATCATAGGCTTCCGTAACAAAACAAAAATCTCCCTGGCCCACTTTTTACCCCAAGGAACACCCACACCCCTAATCCCCATGCTAGTAATCATTGAGACAATCAGCCTATTTATCCAACCGATAGCACTAGCCGTACGACTAACGGCAAACATCACGGCAGGACACCTACTAATGCACCTAATCGGAGGAGCAACCCTCGCATTAATAAACATCAGCATAACCACCGCCCTTATCACGTTCATCATCCTAGTCTTACTAACAGCTCTAGAGTTTGCCGTTGCCATAATCCAAGCGTACGTCTTCACCCTACTAGTAAGTCTATACTTACACGATAACACAT

hippo_orf1 = nt2aa(hippo_orf1)

hippo_orf1 =

MNENLFASFITPTILGLPLVTLIIMFPSILFPAPTRLITNRLVSIQQ*LIQLVSKQIMNIHNHKGQT*TLILISLILFIGSTNLLGLLPHSFTPTTQLSINLGIAIPL*AGTVIIGFRNKTKISLAHFLPQGTPTPLIPMLVIIETISLFIQPIALAVRLTANITAGHLLMHLIGGATLALINISITTALITFIILVLLTALEFAVAIIQAYVFTLLVSLYLHDNT

swalign(hippo_orf1, hippo_orf1)

ans =

  453.6667

[score, localAlignment] = swalign(hippo_orf1, hippo_orf1)

score =

  453.6667


localAlignment =

MNENLFASFITPTILGLPLVTLIIMFPSILFPAPTRLITNRLVSIQQ*LIQLVSKQIMNIHNHKGQT*TLILISLILFIGSTNLLGLLPHSFTPTTQLSINLGIAIPL*AGTVIIGFRNKTKISLAHFLPQGTPTPLIPMLVIIETISLFIQPIALAVRLTANITAGHLLMHLIGGATLALINISITTALITFIILVLLTALEFAVAIIQAYVFTLLVSLYLHDNT
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MNENLFASFITPTILGLPLVTLIIMFPSILFPAPTRLITNRLVSIQQ*LIQLVSKQIMNIHNHKGQT*TLILISLILFIGSTNLLGLLPHSFTPTTQLSINLGIAIPL*AGTVIIGFRNKTKISLAHFLPQGTPTPLIPMLVIIETISLFIQPIALAVRLTANITAGHLLMHLIGGATLALINISITTALITFIILVLLTALEFAVAIIQAYVFTLLVSLYLHDNT

showalignment(localAlignment)
hippo_orf2 = hippo(orfs_valid(1).Start(2) : orfs_valid(1).Stop(2))

hippo_orf2 =

ATGTTAAAATATATTATCCCAACTATCATACTAATACCTCTGACCTGAATATCAAAAAATAGCATAATCTGAACCAACACTACAGCCCACAGTCTGTTAATCAGCTTCACAAGCCTACTCCTCCTAAACCAATTCAACGACAATAGCCTAAACTTCTCACCAATGTTCTTCTCTGACCCCCTATCTACCCCCCTCCTAATCCTAACAATATGGCTCCTACCCCTAATACTAATAGCCAGCCAATCACACCTACTTAAAGAACCCCCAACCCGAAAAAAACTGTTTATCACAATACTAGTTACGCTACAAACATTCCTAATCATAACATTCTCAGCCATAGAACTAATCCTGTTCTATATCCTATTTGAAGCCACACTCATCCCAACCCTCATCATTATCACCCGATGAGGTAACCAAACAGAGCGCCTTAACGCAGGCCTTTACTTTCTATTCTATACCCTAATAGGATCTCTCCCCCTCCTAGTAGCACTGATTTATATCCAAAATATCACAGGATCCCTAAACTTCCTAATACTCCAATATTGAACCCAAGCCGTATCCAACTCTTGATCCAACGTCTTCTTATGACTAGCATGTATAATAGCTTTCATGGTAAAAATACCCCTCTACGGCCTCCACCTATGACTACCTAAAGCACATGTAGAAGCACCCATCGCCGGCTCAATAGTCCTAGCCGCCATTCTACTAAAACTAGGAGGGTATGGCATACTACGTATCACAACTATCCTAAACCCCCTAACAGAAATAATGGCATATCCGTTCATCATACTCTCCCTATGAGGGATAATCATGACTAGCTCCATCTGCCTGCGCCAAACAGACCTGAAATCACTCATCGCATACTCTTCCGTAAGCCACATAGCACTCGTCATTGTAGCAATCCTCATCCAGACCCCATGAAGCTACATAGGAGCAACAGCCCTGATAATCGCCCATGGCCTCACATCATCCATACTGTTCTGCCTAGCAAATTCAAACTACGAGCGAATTCACAGCCGAACAATAATCTTGGCTCGAGGACTACAAACACTCCTTCCACTAATAGCCGCCTGATGACTACTAGCAAGTCTGACAAACCTGGCTCTACCACCCTCTATCAACCTCGTCGGAGAACTACTAGTAATCATGTCCTCTTTCTCATGATCAAATATTACCATTATCCTAATAGGAACCAACATCATAATCACCGCCCTATACTCACTATACATACTAACCACCACTCAACGCGGCAAGTACACCCATCACATCAACAACATTACACCCTCATTTACACGAGAAAACGCCCTAATAGCACTCCACATCCTACCCCTCCTACTACTATCCCTAAATCCCAAAATCATCCTAGGCCCCCTCTACTGTAAGCATAGTTTAA

hippo_orf2 = nt2aa(hippo_orf2)

hippo_orf2 =

MLKYIIPTIILIPLT*ISKNSII*TNTTAHSLLISFTSLLLLNQFNDNSLNFSPMFFSDPLSTPLLILTIWLLPLILIASQSHLLKEPPTRKKLFITILVTLQTFLIITFSAIELILFYILFEATLIPTLIIITR*GNQTERLNAGLYFLFYTLIGSLPLLVALIYIQNITGSLNFLILQY*TQAVSNS*SNVFL*LACIIAFMVKIPLYGLHL*LPKAHVEAPIAGSIVLAAILLKLGGYGILRITTILNPLTEIMAYPFIILSL*GIIMTSSICLRQTDLKSLIAYSSVSHIALVIVAILIQTP*SYIGATALIIAHGLTSSILFCLANSNYERIHSRTIILARGLQTLLPLIAA**LLASLTNLALPPSINLVGELLVIMSSFS*SNITIILIGTNIIITALYSLYILTTTQRGKYTHHINNITPSFTRENALIALHILPLLLLSLNPKIILGPLYCKHSL

hippo_porf1 = nt2aa(hippo(orfs_valid(1).Start(1) : orfs_valid(1).Stop(1)))

hippo_porf1 =

MNENLFASFITPTILGLPLVTLIIMFPSILFPAPTRLITNRLVSIQQ*LIQLVSKQIMNIHNHKGQT*TLILISLILFIGSTNLLGLLPHSFTPTTQLSINLGIAIPL*AGTVIIGFRNKTKISLAHFLPQGTPTPLIPMLVIIETISLFIQPIALAVRLTANITAGHLLMHLIGGATLALINISITTALITFIILVLLTALEFAVAIIQAYVFTLLVSLYLHDNT

hippo_porf2 = nt2aa(hippo(orfs_valid(1).Start(2) : orfs_valid(1).Stop(2)))

hippo_porf2 =

MLKYIIPTIILIPLT*ISKNSII*TNTTAHSLLISFTSLLLLNQFNDNSLNFSPMFFSDPLSTPLLILTIWLLPLILIASQSHLLKEPPTRKKLFITILVTLQTFLIITFSAIELILFYILFEATLIPTLIIITR*GNQTERLNAGLYFLFYTLIGSLPLLVALIYIQNITGSLNFLILQY*TQAVSNS*SNVFL*LACIIAFMVKIPLYGLHL*LPKAHVEAPIAGSIVLAAILLKLGGYGILRITTILNPLTEIMAYPFIILSL*GIIMTSSICLRQTDLKSLIAYSSVSHIALVIVAILIQTP*SYIGATALIIAHGLTSSILFCLANSNYERIHSRTIILARGLQTLLPLIAA**LLASLTNLALPPSINLVGELLVIMSSFS*SNITIILIGTNIIITALYSLYILTTTQRGKYTHHINNITPSFTRENALIALHILPLLLLSLNPKIILGPLYCKHSL

[score, localAlignment] = swalign(hippo_orf2, hippo_orf2)

score =

  923.3333


localAlignment =

MLKYIIPTIILIPLT*ISKNSII*TNTTAHSLLISFTSLLLLNQFNDNSLNFSPMFFSDPLSTPLLILTIWLLPLILIASQSHLLKEPPTRKKLFITILVTLQTFLIITFSAIELILFYILFEATLIPTLIIITR*GNQTERLNAGLYFLFYTLIGSLPLLVALIYIQNITGSLNFLILQY*TQAVSNS*SNVFL*LACIIAFMVKIPLYGLHL*LPKAHVEAPIAGSIVLAAILLKLGGYGILRITTILNPLTEIMAYPFIILSL*GIIMTSSICLRQTDLKSLIAYSSVSHIALVIVAILIQTP*SYIGATALIIAHGLTSSILFCLANSNYERIHSRTIILARGLQTLLPLIAA**LLASLTNLALPPSINLVGELLVIMSSFS*SNITIILIGTNIIITALYSLYILTTTQRGKYTHHINNITPSFTRENALIALHILPLLLLSLNPKIILGPLYCKHSL
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MLKYIIPTIILIPLT*ISKNSII*TNTTAHSLLISFTSLLLLNQFNDNSLNFSPMFFSDPLSTPLLILTIWLLPLILIASQSHLLKEPPTRKKLFITILVTLQTFLIITFSAIELILFYILFEATLIPTLIIITR*GNQTERLNAGLYFLFYTLIGSLPLLVALIYIQNITGSLNFLILQY*TQAVSNS*SNVFL*LACIIAFMVKIPLYGLHL*LPKAHVEAPIAGSIVLAAILLKLGGYGILRITTILNPLTEIMAYPFIILSL*GIIMTSSICLRQTDLKSLIAYSSVSHIALVIVAILIQTP*SYIGATALIIAHGLTSSILFCLANSNYERIHSRTIILARGLQTLLPLIAA**LLASLTNLALPPSINLVGELLVIMSSFS*SNITIILIGTNIIITALYSLYILTTTQRGKYTHHINNITPSFTRENALIALHILPLLLLSLNPKIILGPLYCKHSL

hippo_porf3 = nt2aa(hippo(orfs_valid(1).Start(3) : orfs_valid(1).Stop(3)))

hippo_porf3 =

MEFSIWYIHSDPHINQFFKYLLLFLITIIVLVTANNLFQLFIG*EGVGIMSFLLIG**HGRTDANTAAIQAILYNRIGDVGFIIAIA*FLSNLNT*DIQQIFIINPTHSNLPLIGLILAATGKSAQFGLHP*LPSAMEGPTPVSALLHSSTIVVAGVFLLIRFYPMIENNNLIQTITICLGAITTLFTAICALTQNDIKKIIAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFICSGSIIHNLNNEQDIRKIGGLFKTIPFTTTTLIVGSIALTGVPFLTGFYSKDLIIEAANTSYSNA*ALLITLVATSLTAVYSTRIIFFALLGHPRFPTSTLINENNPLLLNSLKRLIAGSIFAGFILSHNLPPITTPLITIPPYLKMTALAVTILGFTLAFEITLNTQNLKHKHPTNSFKFSTLLGYFPTIIHRLPPHLSLTASQKLASSLLDSA*LENILPKSIAHAQLKLSTLVSNQKGLIKIYFLSFLITIPLSLILFNPHA

[score, localAlignment] = swalign(hippo_orf3, hippo_orf3)
Undefined function or variable 'hippo_orf3'.
 
Did you mean:
[score, localAlignment] = swalign(hippo_porf3, hippo_porf3)

score =

   1.0487e+03


localAlignment =

MEFSIWYIHSDPHINQFFKYLLLFLITIIVLVTANNLFQLFIG*EGVGIMSFLLIG**HGRTDANTAAIQAILYNRIGDVGFIIAIA*FLSNLNT*DIQQIFIINPTHSNLPLIGLILAATGKSAQFGLHP*LPSAMEGPTPVSALLHSSTIVVAGVFLLIRFYPMIENNNLIQTITICLGAITTLFTAICALTQNDIKKIIAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFICSGSIIHNLNNEQDIRKIGGLFKTIPFTTTTLIVGSIALTGVPFLTGFYSKDLIIEAANTSYSNA*ALLITLVATSLTAVYSTRIIFFALLGHPRFPTSTLINENNPLLLNSLKRLIAGSIFAGFILSHNLPPITTPLITIPPYLKMTALAVTILGFTLAFEITLNTQNLKHKHPTNSFKFSTLLGYFPTIIHRLPPHLSLTASQKLASSLLDSA*LENILPKSIAHAQLKLSTLVSNQKGLIKIYFLSFLITIPLSLILFNPHA
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MEFSIWYIHSDPHINQFFKYLLLFLITIIVLVTANNLFQLFIG*EGVGIMSFLLIG**HGRTDANTAAIQAILYNRIGDVGFIIAIA*FLSNLNT*DIQQIFIINPTHSNLPLIGLILAATGKSAQFGLHP*LPSAMEGPTPVSALLHSSTIVVAGVFLLIRFYPMIENNNLIQTITICLGAITTLFTAICALTQNDIKKIIAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFICSGSIIHNLNNEQDIRKIGGLFKTIPFTTTTLIVGSIALTGVPFLTGFYSKDLIIEAANTSYSNA*ALLITLVATSLTAVYSTRIIFFALLGHPRFPTSTLINENNPLLLNSLKRLIAGSIFAGFILSHNLPPITTPLITIPPYLKMTALAVTILGFTLAFEITLNTQNLKHKHPTNSFKFSTLLGYFPTIIHRLPPHLSLTASQKLASSLLDSA*LENILPKSIAHAQLKLSTLVSNQKGLIKIYFLSFLITIPLSLILFNPHA

hippo_porf4 = nt2aa(hippo(orfs_valid(2).Start(1) : orfs_valid(2).Stop(1)))

hippo_porf4 =

MFIINTLILVAPILLAIAFLTLVERKILGYIQLRKGPNIVGPYGLLQPFADAIKLFTKEPLRPSTSSVSIFIIAPILALTLALTI*IPLPIPYPLININLGVLFILAISSLAVYSIL*SG*ASNSKYALIGALRAVAQTISYEVTLAIILLSILLINGSFTLSTLITTQEKL*LIFPS*PLAII*FISTLAETNRAPFDLTEGESELVSGFNVEYAAGPFAIFFIAEYINIIIINAFTTVLFLGAYHNPYLPELYTINFTIKTLLLTISFL*IRASYPRFRYDQLMHLL*KSFLPLTLALCI*HVSLPIITSSIPPQT

[score, localAlignment] = swalign(hippo_porf4, hippo_porf4)

score =

  635.3333


localAlignment =

MFIINTLILVAPILLAIAFLTLVERKILGYIQLRKGPNIVGPYGLLQPFADAIKLFTKEPLRPSTSSVSIFIIAPILALTLALTI*IPLPIPYPLININLGVLFILAISSLAVYSIL*SG*ASNSKYALIGALRAVAQTISYEVTLAIILLSILLINGSFTLSTLITTQEKL*LIFPS*PLAII*FISTLAETNRAPFDLTEGESELVSGFNVEYAAGPFAIFFIAEYINIIIINAFTTVLFLGAYHNPYLPELYTINFTIKTLLLTISFL*IRASYPRFRYDQLMHLL*KSFLPLTLALCI*HVSLPIITSSIPPQT
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MFIINTLILVAPILLAIAFLTLVERKILGYIQLRKGPNIVGPYGLLQPFADAIKLFTKEPLRPSTSSVSIFIIAPILALTLALTI*IPLPIPYPLININLGVLFILAISSLAVYSIL*SG*ASNSKYALIGALRAVAQTISYEVTLAIILLSILLINGSFTLSTLITTQEKL*LIFPS*PLAII*FISTLAETNRAPFDLTEGESELVSGFNVEYAAGPFAIFFIAEYINIIIINAFTTVLFLGAYHNPYLPELYTINFTIKTLLLTISFL*IRASYPRFRYDQLMHLL*KSFLPLTLALCI*HVSLPIITSSIPPQT

hippo_porf5 = nt2aa(hippo(orfs_valid(2).Start(2) : orfs_valid(2).Stop(2)))

hippo_porf5 =

MIVITSSH*LLT*TGFEINMLAIIPIIMKSPNPRATEASVKYFITQATASMLLILAVIINLLYSGQ*TVIKILNPTASIIITIALAIKLGLSPFHF*VPEVTQGIPLTAGLILLT*QKLAPLSILYQISPSINPNLILTISMLSILVGG*GGLNQTQLRKIIAYSSIAHMGWIAAILIYNPTITILNLTIYLITTFTMFTIFALNSTTTTLSLSHT*NKTPIITTLILTILLSIGGLPPLTGFVPK*IIIQEITKNDSIILPTLIAIIALPNLYFYIRLTYSTALTIFPSSNNIKIK*QFEASKHKTLLPTIIILSTILLPLTPILVVLD

[score, localAlignment] = swalign(hippo_porf5, hippo_porf5)

score =

  648.3333


localAlignment =

MIVITSSH*LLT*TGFEINMLAIIPIIMKSPNPRATEASVKYFITQATASMLLILAVIINLLYSGQ*TVIKILNPTASIIITIALAIKLGLSPFHF*VPEVTQGIPLTAGLILLT*QKLAPLSILYQISPSINPNLILTISMLSILVGG*GGLNQTQLRKIIAYSSIAHMGWIAAILIYNPTITILNLTIYLITTFTMFTIFALNSTTTTLSLSHT*NKTPIITTLILTILLSIGGLPPLTGFVPK*IIIQEITKNDSIILPTLIAIIALPNLYFYIRLTYSTALTIFPSSNNIKIK*QFEASKHKTLLPTIIILSTILLPLTPILVVLD
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MIVITSSH*LLT*TGFEINMLAIIPIIMKSPNPRATEASVKYFITQATASMLLILAVIINLLYSGQ*TVIKILNPTASIIITIALAIKLGLSPFHF*VPEVTQGIPLTAGLILLT*QKLAPLSILYQISPSINPNLILTISMLSILVGG*GGLNQTQLRKIIAYSSIAHMGWIAAILIYNPTITILNLTIYLITTFTMFTIFALNSTTTTLSLSHT*NKTPIITTLILTILLSIGGLPPLTGFVPK*IIIQEITKNDSIILPTLIAIIALPNLYFYIRLTYSTALTIFPSSNNIKIK*QFEASKHKTLLPTIIILSTILLPLTPILVVLD

clear hippo_orf1
clear hippo_orf2
human_porf1 = nt2aa(human(orfs_human(1).Start(1) : orfs_human(1).Stop(1)))

human_porf1 =

MANLLLLIVPILIAMAFLMLTERKILGYIQLRKGPNVVGPYGLLQPFADAIKLFTKEPLKPATSTITLYITAPTLALTIALLL*TPLPIPNPLVNLNLGLLFILATSSLAVYSIL*SG*ASNSNYALIGALRAVAQTISYEVTLAIILLSTLLISGSFNLSTLITTQEHL*LLLPS*PLAII*FISTLAETNRTPFDLAEGESELVSGFNIEYAAGPFALFFIAEYTNIIIINTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFL*IRTAYPRFRYDQLIHLL*KNFLPLTLALLI*YVSIPITISSIPPQT

[score, localAlignment] = swalign(hippo_porf1, human_porf1)

score =

   22.3333


localAlignment =

ILISLILFIGSTNLLGLLPHS-FTPTTQLSINLGIAIPL-*AGTVIIGFRNKTKISLAHFLPQGTPTPLIPM-LVIIETISLFIQPI-ALAV-RLTANITAGHLLMHLIGG-ATLA-LINISITTALITFIILVLLTALEFAVAI-IQAYVFTLLVS
||  : |  |  |::|  |:: : | :: :|:|    ||  | ::|  : :   ::|:  |   || | ||  || ::   |||    :|||  :  :  |::  : |||:  ::|  |:  :| |:| :  |::  :::::: |  | ::  || |
ILGYIQLRKG-PNVVG--PYGLLQPFAD-AIKLFTKEPLKPATSTITLYITAPTLALTIALLL*TPLP-IPNPLVNLNLGLLFILATSSLAVYSIL*SG*ASNSNYALIGALRAVAQTISYEVTLAIILLSTLLISGSFNLSTLITTQEHL*LLLPS

funcs_mt1
funcs_mt1
funcs_mt1
funcs_mt1
[score, localAlignment] = nwalign(hippo_porf1, human_porf1)

score =

 -102.6667


localAlignment =

M-N------ENLFA-SF--ITP-TILG---L---P-LV-TL-IIM-F-PSI-LF---P-AP-TRLITNRLVSIQQ*L-IQ-LV-SK-QIMN-IHN-HKGQT*TLILISLILF-I---G-STNLLGLLPHSFTPTTQ-LSINLGIAIPL*AGTVIIG-FRNKTKISLAHFLPQGTPT-PLIPMLVIIETISLFIQ-PIALA---VRL-TA-NI--TAGHLLMHLIGGAT-LALIN-I--SI---TT--AL---I--T-FI---ILV--L-L---TAL-EF---A-VAII-QAYV-FT--LL---VSL-Y-LHD---NT
| |        |:| :|  :|   |||   |   | :|    ::: |  :| ||   |  | |  ||  :::    | |  |: :   | | : | : |    |   || :: |   | ::|    |  ::  ::| :| :: :|| | :  :| | |  :| |:  : |    |: ||  :: :| |::   : |: ||    :| :: ||  :|| : : :|:  | : :|| :  :|   ||  ||   :  | |:   :|:  | |   ||  :|     : :: : :: :|  ||   ||:   : :   :|
MANLLLLIVPILIAMAFLMLTERKILGYIQLRKGPNVVGPYGLLQPFADAIKLFTKEPLKPATSTITLYITAPTLALTIALLL*TPLPIPNPLVNLNLGLLFILATSSLAVYSIL*SG*ASNSNYALIGALRAVAQTISYEVTLAIILLSTLLISGSFNLSTLITTQEHL*LLLPS*PL-AII*FISTLAETNRTPFDLAEGESELVSGFNIEYAAGPFALFFIAEYTNIIIINTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFL*IRTAYPRFRYDQLIHLL*KNFLPLTLALLI*YVSIPITISSIPPQT

n = 1000;
globalscores = zeros(n,1);
hippoLen = length(hippo_porf1);
for i = 1:n
    permutedSequence = hippo_porf1(randperm(hippoLen));
    globalscores(i) = nwalign(human_porf1,permutedSequence);
end

figure
buckets = ceil(n/20);
hist(globalscores,buckets)
hold on;
stem(score,1,'c')
title('Determining Alignment Significance: Hippo vs Human Proteins');
xlabel('Score'); ylabel('Number of Sequences');


parmhat = evfit(globalscores)

parmhat =

 -106.4099    8.1155
 
 p = 1 - evcdf(score,parmhat(1),parmhat(2));
 
 human_porf2 = nt2aa(human(orfs_valid(1).Start(2) : orfs_valid(1).Stop(2)))

human_porf2 =

TLYPPPASLSP*NSS**LLPSYYLI*KLPSFYPYHEPYKQLTCH**LCHPSY*SSS*P*VWPMSDYKKD*TEPNWYIV*TKRMISTH*IMIIIFTKCPSFT*ILY*HLPSHF*EY*YIAHTSCPPYYA*KE*YYRCSL*LLS*PSTPTPS*PILCLLPY*SLPPAKQRWA*PY*SQSPTHMA*TTYIT*TYSNAKTNRPNNYITTTDMTFQKTHNLNQHNHPQPNY*HHPSTIF*PNQQQPI*LFPNLFLRPPNNPPPNTNYLTPTPHNHGKPTPLIQ*TTITKKTLPLYTNLPTNLLNYNIHSHRTNHILYLLRNHTYPHLGYHHPMRQPARTPERRHILPILHPSRLPSPTHRTNLHSQHPRLTKHSTTHPHCPRTIKLLSQQLNMTSLHNSFYSKDTSLRTPLMTP*SPCRSPHRWVNSTCRSTLETRRLWHNTPHTHSQPPDKTHSLPLPCTIPMRHNYN

[score, globalAlignment] = nwalign(hippo_porf2, human_porf2)

score =

  -87.6667


globalAlignment =

MLKYIIPTIILIP--LT-*-ISKNSII*TNTTAHSLLISFTSLLLLNQFNDNSL-NFS-P-MF-FSDPLSTPLLILTIWLLPLILIASQSH-LLKEPPTRKKLFITILV-TLQT-FLIITFSA-IELILFYILFEATLIPTLIIITR*GNQTERLNAGLYFLFYTLIGSLPLLVALIY-IQNITG-S-LNFLILQY*TQAVSNS*SNVFL*LACIIAFMVKIPLYGLH-L*LPK--AH---VEAPIAGS-I-VLAAILLKLGGYGILRITTILNPLTEIMAYPFIILSL*GIIMTSSICLRQTDLKS-LIAYSSVSHIALVIVAILIQTP*SYIGATALIIAHGLTSSILFCLANSNYERIHSRTIILARGLQTLLPLIAA**LLASLTNLALPPSINLVGELLVIMSSFS*SNITIILIGTNIIITALYSLYILTTTQRGKYTHHINNITPSFTRENALIALHILPLLLLSLNP-KIILGPLYC----KHSL-
 | |  |   | |   : | : :  :||   : :     : :|     :   |  : | | :: :||  :        | :       ::| ::    |:   |  ||   | : |    : |      :|   |     :| :::  :: |      | :| |  :       |  |  |: |  :  :::   | ::| :|   | ::  :  ::|: |    : |    |:   |   :  |   : | ::  ::|:  :      |: |:|  :  : |  ::: |  |  ::: |  |:| : |: |:  || :  |: :| :    ::|    :   : |      |   :  |: | |     :|::  | ::      | |:   | :|:|::: |  |:|:   | ::    |: : | |::    :   |: : : :|:   | | |:  :  |  |    |  | |    || |    :|:  
TL-Y-PPPASLSP*NSS**LLPSYYLI*KLPSFYPYHEPYKQLTCH**LCHPSY*SSS*P*VWPMSD-YKKD-*TEPNWYIV*TKRMISTH*IMIIIFTKCPSFT*ILY*HLPSHF*EY*YIAHTSCPPYYA*KE*YYRCSL*LLS*PSTPTPS*PI-LCLLPY*SLPPAKQRWA*PY*SQSPTHMA*TTYIT*TY-SNAKTNR-PNNYI-TTTDMTFQ-KTHNLNQHNHPQPNY*HHPSTIF*PNQQQPI*LFPNLFLRPPNNPPPN-TNYLTPTPHNHGKPTPLIQ-*TTITKKTLPL-YTNLPTNLLNYNIHSHRTNHILYLLRNHTYPHLGYHHPMRQPARTPERRHILPILHPSRLPSPT--HRTNLHSQHPRLTK----HSTTHPHCPRTIKLLSQQL-NMTSL--HN-SFYSKDTS-LRTPLMT--P*SPC-RSPH-RWVNSTCRS-TLETRRL-WHNTP-HTHSQPPDKTHSLPLPCTIPMRHNYN

funcs_mt1
funcs_mt1
funcs_mt1

score =

  -403


globalAlignment =

MLKYIIPTIILIPLT*ISKNSII*TNTTAHSLLISFTSLLLLNQFNDNSLNFSPMFFSDPLSTPLLILTIWLLPLILIASQSHLLKEPPTRKKLFITILVTLQTFLIITFSAIELILFYILFEATLIPTLIIITR*GNQTERLNAGLYFLFYTLIGSLPLLVALIYIQNITGSLNFLILQY*TQAVSNS*SNVFL*LACIIAFMVKIPLYGLHL*LPKAHVEAPIAGSIVLAAILLKLGGYGILRITTILNPLTEIMAYPFIILSL*GIIMTSSICLRQTDLKSLIAYSSVSHIALVIVAILIQTP*SYIGATALIIAHGLTSSILFCLANSNYERIHSRTIILARGLQTLLPLIAA**LLASLTNLALPPSINLVGELLVIMSSFS*SNITIILIGTNIIITALYSLYILTTTQRGKYTHHINNITPSFTRENALIALHILPLLLLSLNPKIILGPLYCKHSL
| : :: ::|  | |      |:  :  | :  :     |::       |      |  |   | | | |   |    :|: :|: :  :|  | ||   | |  | | ::: :  :  |    | | :    |: |    | : :|  :   |: |  |:: :|   |: | |  |  :: | ::: |  |    |    | |: ||     ||          ||:::   |     |:   | : | | |    |: |  |   |   | :  |    :  |:    :  |  :| ::||| |   |   || : : ||       |  |   |   |   | | : || :     |::| |:||:  :|||    |:|   |  : | :||: | |:  |: :| |    | : |      :   ||    |  | :      ||        |:  :: 
MNENLFASFI-AP-T------IL--GLPA-A--V-----LII-------L------F--P---P-L-L-I---P----TSK-YLI-N--NR--L-IT---T-QQ*L-IKLTS-KQ-M--I----T-IHN----TK-G----R-T*SL--I---LV-S--LII-FIATTNLLG-L--LPHSF-TPTTQLS-IN----L----A-MA-IP-----L*----------AGAVI---I-----GF---R-SKIKNAL----AH-F--LPQ-G---TPTP-L--IPI--LV----I--IE-TI-SLLIQ-P---I---ALAV-R-LT-------A--N---I---T---A-G-H-LL-M--H--LIGS-TTLAI-STINL-PSTLII---F--T-I-LILL-T-ILEIAV-AL-I----Q-A-Y------V---FT----L--L-V------SL-----Y--LH-DNT-


p =

    0.3007

funcs_mt1

score =

  492.3333


globalAlignment =

----------------------------------------------------M----E------------------------FS---I-------W--------YIHSDPHINQFFKYLLLFLITIIVLVTANNLFQLFIG*EGVGIMSFLLIG**HGRTDANTAAIQAILYNRIGDVGFIIAIA*FLSNLNT*DIQQIFIINPTHSNLPLIGLILAATGKSAQFGLHP*LPSAMEGPTPVSALLHSSTIVVAGVFLLIRFYPMIENNNLIQTITICLGAITTLFTAICALTQNDIKKIIAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFICSGSIIHNLNNEQDIRKIGGLFKTIPFTTTTLIVGSIALTGVPFLTGFYSKDLIIEAANTSYSNA*ALLITLVATSLTAVYSTRIIFFALLGHPRFPTSTLINENNPLLLNSLKRLIAGSIFAGFILSHNLPPITTPL-ITIPPYLKMTALAVTILGFTLAFEITLNTQNLKHKHPTNSFKFSTLLGYFPTIIHRLPPHLSLTASQKLASSLLDSA*LENILPKSIAHAQLKLSTLVSNQKGLIKIYFLSFLITIPLSLILFNPHA
                                                    |    |                        ||   |       |        ||:|||:|||||||||:|||||::|||||||||||||||||||:|||||:||::|:|||||||||||||||||:|||:|:|||: : |:|| ||| ::| : |  ||:||:|||:|||||:|||||||||:|||||||||||||||||||:||||||:|: ||: ||||:|:|||||||||:|:|||||||||||:||||||||||||||||||||||||||||||||||||||:|||||||||||||||||||||:||||:|:|:| :||:||:|:|||||||||| |||:|| ||:||||| |||:|||||::||||||:::| |:||||| | |||||| ||| :||| |||:||||::::|: | ::|:  ||| |||:||||||:||:  |::::  |::|| | |  :| ||::||::|:| ||  |:|:| :||:|   ||| :|||::|||:|:: |:: | ::|:|||:||:|||||:: : |:|:|:   :
MHTTITTLTLTSLIPPILTTLVNPNKKNSYPHYVKSIVASTFIISLFPTTIFMCLDQEVIISN*H*ATTQTTQLSLSFKLDYFSIIFIPVALFVTWSIIEFSL*YINSDPNINQFFKYLLIFLITILILVTANNLFQLFIG*EGVGIISFLLIS**YARADANTAAIQAILYNRIGDIGFILALA*FILHSNS*DPQQIALLNANPSLTPLLGLLLAAAGKSAQLGLHP*LPSAIEGPTPVSALLHSSTIVVAGIFLLIRFHPLAENSPLIQTLTLCLGAITTLFAAVCALTQNDIKKIVAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFMCSGSIIHNLNNEQDIRKIGGLLKTIPLTSTSLTIGSLALAGIPFLTGFYSKDHIIETANISYTNA*ALSITLIATSLTSAYSTRIILLTLTGQPRFPTLTNINENNPTLLNPIKRLAAGSLFAGFLITNNISP-ASPFQTTIPLYLKLTALAVTFLGLLTALDLNYLTNKLKIKSPLCTFYFSNILGFYPSITHRTIPYLGLLTSQNLPLLLLDLT*LEKLLPKTISQHQISTSIITSTQKGIIKLYFLSFFFPLILTLLLI---T


p =

     0

funcs_mt1

score =

  492.3333


globalAlignment =

----------------------------------------------------M----E------------------------FS---I-------W--------YIHSDPHINQFFKYLLLFLITIIVLVTANNLFQLFIG*EGVGIMSFLLIG**HGRTDANTAAIQAILYNRIGDVGFIIAIA*FLSNLNT*DIQQIFIINPTHSNLPLIGLILAATGKSAQFGLHP*LPSAMEGPTPVSALLHSSTIVVAGVFLLIRFYPMIENNNLIQTITICLGAITTLFTAICALTQNDIKKIIAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFICSGSIIHNLNNEQDIRKIGGLFKTIPFTTTTLIVGSIALTGVPFLTGFYSKDLIIEAANTSYSNA*ALLITLVATSLTAVYSTRIIFFALLGHPRFPTSTLINENNPLLLNSLKRLIAGSIFAGFILSHNLPPITTPL-ITIPPYLKMTALAVTILGFTLAFEITLNTQNLKHKHPTNSFKFSTLLGYFPTIIHRLPPHLSLTASQKLASSLLDSA*LENILPKSIAHAQLKLSTLVSNQKGLIKIYFLSFLITIPLSLILFNPHA
                                                    |    |                        ||   |       |        ||:|||:|||||||||:|||||::|||||||||||||||||||:|||||:||::|:|||||||||||||||||:|||:|:|||: : |:|| ||| ::| : |  ||:||:|||:|||||:|||||||||:|||||||||||||||||||:||||||:|: ||: ||||:|:|||||||||:|:|||||||||||:||||||||||||||||||||||||||||||||||||||:|||||||||||||||||||||:||||:|:|:| :||:||:|:|||||||||| |||:|| ||:||||| |||:|||||::||||||:::| |:||||| | |||||| ||| :||| |||:||||::::|: | ::|:  ||| |||:||||||:||:  |::::  |::|| | |  :| ||::||::|:| ||  |:|:| :||:|   ||| :|||::|||:|:: |:: | ::|:|||:||:|||||:: : |:|:|:   :
MHTTITTLTLTSLIPPILTTLVNPNKKNSYPHYVKSIVASTFIISLFPTTIFMCLDQEVIISN*H*ATTQTTQLSLSFKLDYFSIIFIPVALFVTWSIIEFSL*YINSDPNINQFFKYLLIFLITILILVTANNLFQLFIG*EGVGIISFLLIS**YARADANTAAIQAILYNRIGDIGFILALA*FILHSNS*DPQQIALLNANPSLTPLLGLLLAAAGKSAQLGLHP*LPSAIEGPTPVSALLHSSTIVVAGIFLLIRFHPLAENSPLIQTLTLCLGAITTLFAAVCALTQNDIKKIVAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFMCSGSIIHNLNNEQDIRKIGGLLKTIPLTSTSLTIGSLALAGIPFLTGFYSKDHIIETANISYTNA*ALSITLIATSLTSAYSTRIILLTLTGQPRFPTLTNINENNPTLLNPIKRLAAGSLFAGFLITNNISP-ASPFQTTIPLYLKLTALAVTFLGLLTALDLNYLTNKLKIKSPLCTFYFSNILGFYPSITHRTIPYLGLLTSQNLPLLLLDLT*LEKLLPKTISQHQISTSIITSTQKGIIKLYFLSFFFPLILTLLLI---T


p =

     0

funcs_mt1

score =

 -129.6667


globalAlignment =

MFIINTLILVAPILLAIAFLTLVERKILGYIQLRKGPNIVGPYGLLQPFADAIKLFTKEPLRPSTSSVSIFIIAPILALTLALTI*IPLPIPYPLININLGVLFILAISSLAVYSIL*SG*ASNSKYALIGALRAVAQTISYEVTLAIILLSILLINGSFTLSTLITTQEKL*LIFPS*PLAII*FISTLAETNRAPFDLTEGESELVSGFNVEYAAGPFAIFFIAEYINIIIINAFTTVLFLGAYHNPYLPELYTINFTIKTLLLTISFL*IRASYPRFRYDQLMHL--L*KSFLPLTLALCI*HVSLPIITSSIPPQT-
|   :    :| : |  |   ::| : |  | ::    |:  : |:  |   : |::   |   |::::   |:   |  :  |:|  : :|   | | | || | |: || :  |  :    |:  :|  :::::::   |  |      : |::| |:     |       |  |  ||    |:   ::  |  :|:   ::::|  : :|  | |: | :|:   :  ::::  :|  ||  ::  :|    :|   :|   |  |    |: |   | |  ::    :||:|::| | |    | |:  : |  | 
M--AH----AAQVGLQDATSPIIE-E-L--ITFHDHALII-IF-LIC-F--LV-LYALF-L-TLTTKLTNTNISD--AQEIE-TV*-TI-LP-A-I-I-L-VL-I-ALPSLRILYI--TD-EVNDP-SL--TIKSIGHQW-Y-*TYEYTDYGGLIFN-SY-----I-------L--P--PL----FLEP-GDL-RL-LDV---DNRVV--LPIE--A-PIRI-IITSQ-D--VLHS-*AVPTLG-LKTDAIP--GRLN---QT---T--F---TATRPGVYYGQCSEICGANHSFMPIVLEL-I---PLKIFEIG-PVFTL


p =

    0.7253

funcs_mt1

score =

 -129.6667


globalAlignment =

MFIINTLILVAPILLAIAFLTLVERKILGYIQLRKGPNIVGPYGLLQPFADAIKLFTKEPLRPSTSSVSIFIIAPILALTLALTI*IPLPIPYPLININLGVLFILAISSLAVYSIL*SG*ASNSKYALIGALRAVAQTISYEVTLAIILLSILLINGSFTLSTLITTQEKL*LIFPS*PLAII*FISTLAETNRAPFDLTEGESELVSGFNVEYAAGPFAIFFIAEYINIIIINAFTTVLFLGAYHNPYLPELYTINFTIKTLLLTISFL*IRASYPRFRYDQLMHL--L*KSFLPLTLALCI*HVSLPIITSSIPPQT-
|   :    :| : |  |   ::| : |  | ::    |:  : |:  |   : |::   |   |::::   |:   |  :  |:|  : :|   | | | || | |: || :  |  :    |:  :|  :::::::   |  |      : |::| |:     |       |  |  ||    |:   ::  |  :|:   ::::|  : :|  | |: | :|:   :  ::::  :|  ||  ::  :|    :|   :|   |  |    |: |   | |  ::    :||:|::| | |    | |:  : |  | 
M--AH----AAQVGLQDATSPIIE-E-L--ITFHDHALII-IF-LIC-F--LV-LYALF-L-TLTTKLTNTNISD--AQEIE-TV*-TI-LP-A-I-I-L-VL-I-ALPSLRILYI--TD-EVNDP-SL--TIKSIGHQW-Y-*TYEYTDYGGLIFN-SY-----I-------L--P--PL----FLEP-GDL-RL-LDV---DNRVV--LPIE--A-PIRI-IITSQ-D--VLHS-*AVPTLG-LKTDAIP--GRLN---QT---T--F---TATRPGVYYGQCSEICGANHSFMPIVLEL-I---PLKIFEIG-PVFTL


p =

    0.7354

funcs_mt1

score =

 -177.6667


globalAlignment =

MI-VITSS-H*L-LT*-T--GFE-IN-----M-LAIIPII-M-K-SPN-----PR-ATEA-SVKYFI-TQ---A-T--ASM--L--LILA---VIINLLYSGQ*TVIKILNPTASII--I----T-I-ALA-I-KLGLSPFHF*V-PE-V--TQ-G-IP----L--T--A-G-L-I--L-LT*QKLA-P-L-S-I-L-YQIS-----P--SIN---PNL-I-LTI--SM-L-SI-L-VGG*GGLNQT---Q-LRKIIAYS-SIAHM-G-WI-AAI---------LIYNPTIT-I-LNLTIYLI-T--TFT--MFTIFALNSTTTTL-SLSHT-*NKT-P--II-T----TLI-L-TI--LL-SI-G-GLPPLTGFVPK-*IIIQEIT-KNDSIILPTL-IAIIALPNLYFYIRLTY-S-T-AL-TIFPS-S-NNIKIK*QFEASKHKTLLPTIII-LSTILLPLTPIL-VVLD-
|: :|: :   | ||| :   :  ||     : ::|||:: : : : |     |  :::  ::  :| |      |  ||:  |    |:   : :::| | | ::|  :: |  ||  |    | | :|| | : | :| :: :    :  |  | :|    |  |  : | | |  | || |:|:     : | | | |:     |  :::   |:  :   |  |: | :: | :|| | :  |   : | | |||   :  : |  | ::|         ||   :|: | | :|  || |  :||  :: |:| : |:: |  |:::  ::|    || :    ||: | ::  || |: : :|||  ::: :  :::  :: :| :::|  | | : || :||::    : | |  : :| || : :|  :  ::      :| | ||  :|:    ||       | 
MLKLIVPTIILLPLT*LSKKHII*INTTTHSLIISIIPLLFFNQINNNLFSCSPTFSSDPLTTPLLILTT*LLPLTIMASQRHLSSEPLSRKKLYLSILISLQISLIITFTATELIIFYIFFETTLIPTLAIITR*GNQPERLNAGTYFLFYTLVGSLPLLIALIYTHNTLGSLNILLLTLTAQELSNS*ANNLI*LAYTIAFIVKIPLYGLHL*LPKAHVEAPIAGSIVLAAVLLKLGGYGIIRLTLILNPLTKHIAYPFLVLSL*GIIITSSICLRQTDLKSLIAYSSISHIALVVTAILIQTP*SFTGAVILIIAHGLTSSLLFCLANSNYERTHSRIIILSQGLQTLLPLIAF**LLASLANLALPPTINLLGELSVLVTTFS*SNITLLLTGLNILVTALYSLYIFTTTQWGSLTHHINNIKPSFTRENTLMFIHLSPILLLSLNPDIITGFSSCKYSLTKTSDCESDN


p =

    0.2010

funcs_mt1

score =

 -181.3333


globalAlignment =

MFINR*LFSTNHKDIGTLYLLFGA*AGIAGTGLSLLIRAELGQPGTLLGDDQIYNVVVTAHAFVIIFFIVIPIMIGGFGN*LVPLIIGAPDMAFPRINNISF*LLPPSFLLLLASSMVEAGAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAISQYQTPLFV*SVLITAVLLLLSLPVLAAGITMLLTDRNLNTTFFDPAGGGDPVLYQHLF*FFGHPEVYILILPGFGIISHIVTYYSGKKEPFGYIGIV*AIISIGFLGFIV*AHHIFTVGIDVDTRAYFTSATIIIAIPTGVKVFS*LATLHGGNIKWSPAMM*ALGFIFLFTVGGLTGIVLANSSLDIVLHDTYYVVAHFHYVLSIGAVFAIIGGFVH*FPLFSGYTLNDT*AKIHFVIMFVGVNLTFFPQHFLGLSGMPRRYSDYPDAYTT*NTISSIGSFISLTAVVLIVFII*EAFVSKREVLAVDLTTTNLE*LNGCPPPYHTFEEPAYVNLTSQNKRG
|   |   : |   |  | :  :    :  |  |  | |   : |:|||   | : ::|:  |: : :   |    :|::  :   | : |: :  |  | : |   :  :::   :     | |  :|    |:: ::  : :: |: | |  :::|  |  ||   : :  |  :|  :  |  :::||   || ::| :  | | |:  : :    :: :   :|:|  ::|  | |  :  :|: ::: :  :  : :|:::|:|  :     |::     |  |  ::   ::: | : | |   ::: | :  : : | |    :||  |   |   |:| :|:  | : |  | :|  ||       | :  :|||  ::|||  ::| :   ::  :    ::  :  : ::::| ::    :|  |:|   | ||  : |   |::::| : :| ::||  :|  : :|    | :: :   |:  : |  |       : :     |   
MTPIR---KIN-PLI-KL-INHSL-IDLP-TP-S-NISA-**NFGSLLGACLILQ-ITTG-LFLAMHY-S-PDASTAFSS--IAH-I-TRDVNYG*I--IRY-LHANGASIFFI-CLF-LHIGRG--LY--Y-GSFLYS-ET*NIGII-L-L--LATI--ATAFI-GYV-L--P-*GQI-S--F*GATVITN--LLSAIPYI--G-TDLV--Q*I*GG-YS-V--DSPTL-TRFF-TF-H-FILPFIIAALAALHLLFLHETGSNNPLGITSHS-DKITFHPYYTIKDALGLLLFLLSLMTLTLF-SPD-LLGDPDNYTLANPLNT--PPHIK--PE--*-Y-FLFAYTI--LRS-V-PN-KLGGVL-------A-L--LLSI-LILAII-PILH-ISKQQSIIFRPL-SQ-SLY*L-LAADLLIL-T-*IG--GQP--VS-YP--F-T--IIGQVASVLYFT-TILI--LI-PT-IS----L-IE-NKI-LK-WT-C--P--C----S-I-----N---


p =

    0.7872

funcs_mt1

score =

 -181.3333


globalAlignment =

MFINR*LFSTNHKDIGTLYLLFGA*AGIAGTGLSLLIRAELGQPGTLLGDDQIYNVVVTAHAFVIIFFIVIPIMIGGFGN*LVPLIIGAPDMAFPRINNISF*LLPPSFLLLLASSMVEAGAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAISQYQTPLFV*SVLITAVLLLLSLPVLAAGITMLLTDRNLNTTFFDPAGGGDPVLYQHLF*FFGHPEVYILILPGFGIISHIVTYYSGKKEPFGYIGIV*AIISIGFLGFIV*AHHIFTVGIDVDTRAYFTSATIIIAIPTGVKVFS*LATLHGGNIKWSPAMM*ALGFIFLFTVGGLTGIVLANSSLDIVLHDTYYVVAHFHYVLSIGAVFAIIGGFVH*FPLFSGYTLNDT*AKIHFVIMFVGVNLTFFPQHFLGLSGMPRRYSDYPDAYTT*NTISSIGSFISLTAVVLIVFII*EAFVSKREVLAVDLTTTNLE*LNGCPPPYHTFEEPAYVNLTSQNKRG
|   |   : |   |  | :  :    :  |  |  | |   : |:|||   | : ::|:  |: : :   |    :|::  :   | : |: :  |  | : |   :  :::   :     | |  :|    |:: ::  : :: |: | |  :::|  |  ||   : :  |  :|  :  |  :::||   || ::| :  | | |:  : :    :: :   :|:|  ::|  | |  :  :|: ::: :  :  : :|:::|:|  :     |::     |  |  ::   ::: | : | |   ::: | :  : : | |    :||  |   |   |:| :|:  | : |  | :|  ||       | :  :|||  ::|||  ::| :   ::  :    ::  :  : ::::| ::    :|  |:|   | ||  : |   |::::| : :| ::||  :|  : :|    | :: :   |:  : |  |       : :     |   
MTPIR---KIN-PLI-KL-INHSL-IDLP-TP-S-NISA-**NFGSLLGACLILQ-ITTG-LFLAMHY-S-PDASTAFSS--IAH-I-TRDVNYG*I--IRY-LHANGASIFFI-CLF-LHIGRG--LY--Y-GSFLYS-ET*NIGII-L-L--LATI--ATAFI-GYV-L--P-*GQI-S--F*GATVITN--LLSAIPYI--G-TDLV--Q*I*GG-YS-V--DSPTL-TRFF-TF-H-FILPFIIAALAALHLLFLHETGSNNPLGITSHS-DKITFHPYYTIKDALGLLLFLLSLMTLTLF-SPD-LLGDPDNYTLANPLNT--PPHIK--PE--*-Y-FLFAYTI--LRS-V-PN-KLGGVL-------A-L--LLSI-LILAII-PILH-ISKQQSIIFRPL-SQ-SLY*L-LAADLLIL-T-*IG--GQP--VS-YP--F-T--IIGQVASVLYFT-TILI--LI-PT-IS----L-IE-NKI-LK-WT-C--P--C----S-I-----N---


p =

    0.7807

funcs_mt1

score =

  -113


globalAlignment =

M-AY-P-L--QL-------G---F-QDAVSPIIEEL-LYFHD----H-TLI-IVFLISSL-VLY-I-ITL-ILTTKL-T-H-T-NT-INAQEVETV*TIL-P-AIIL-I---L-IA-LPSLRILYII-DE---INNPSL-TV---KTMGH-Q*Y*S-YEYT-DYEDLNFDSYIV-PTS-----DLKPGDLR-LL-EVDNRVV-L-P-IDVT-VRI--LISSEDVLHS*AVPSLGLKTDA--IP---G--R-LN---QTTLI-STRPGLF-YGQCSEIC-GSNH---S-FMPIVLELVP-LQTFEK*TASLL
| |: | |  ::       :   |  :|:: ||  : : |::    : |:   :   ||| ::: | | | |   :: : : | :| :::  :  :|  | | :||  |   | :: | :| || ||      :|: :|  :   ::: |  |  :   |: :   ||:  ||:  |:     :|: :    || :: |::: | | |  | : :  |      | : |:     |:::  ||   :    ||      || ||   |:  ::  :|    :|   : |:| :: |:  |  :      :|
MLAFIPVLTKKINPRSTEAAIKYFLTQATASIILLIAILFNNILSGQ*TITNTTNQYSSLIIIMAIAIKLGIAPFHF*VPEVTQGTPLTSGLLLLT*QKLAPISIIYQISPSLNVSLLLTLSILSIIAGS*GGLNQTQLRKILAYSSITHIG*IIAVLPYNPNITILNLTIYIILTTTAFLLLNLNSSTTTLLLSRT*NKLT*LTPLIPSTLLSLGGLPPLTGFLPKWAIIEEFTKNNSLIIPTIIATITLLNLYFYLRLIYSTSITLLPISNNVKIK*QFEHTKPTPFLPTLIALTTLLLPISPFILIIL


p =

    0.4468

funcs_mt1

score =

 -448.3333


globalAlignment =

M--------T-HQ---T-H----A-YHIV-NP-S-----P--*P--LTGA--LSALLITS-GLTI*FH---------F-NSLILLTTG---L----VTNI----L--TIYQ**RDVIRE----S--T-F---QGHHT-P-V-VQ-K--GLRY-GI--VL----F---IIS-E--VLF-F-TGFF-*-----A--FYHS-S-LAP--T---PE--LGGC-*PPT-GINP-----L-----NP---LEV-P---LL-N--T--S--------V-LL-A--S-G-VS-IT*AHH--SL-IE-GNRKQILQA-LFITIALGV-YFT-L--LQASEYH--EA---S------FTI-S-DG-VY--GS-------TFFVATGFH-----G-LHVIIGS-T--F-LIVCF-LRQ-L-KFHFT----S-D-----HHF----G----F-EAAA*Y--*H-------FVD---VV-*LFL----Y--------V---SI-Y-*-*G---------------S
|        | |:   | :    |   :: :  |         |  | |   :  :::|: :::| |          | | |: |  |   :    ::||    |  ::     ::| |    :  | :    |::: | : |:    :|:  |:  :|    |   ||: :  ::  : | :| |     |  :  |   ||   |    :  |:     |: | :|     |     :|   : : |   :: :  |  |        : :: |  | | :: |:||||  :: |:  :|  : :| ::|:|  ||  |: |  |::|:::   |   :      ||: :  | |   :|       |::|:: ||     | : :|||:    | |:  : | |   |:|||    : :     :||    |    : :    |  |:       |::   |:  :|:    :        |   ||   |  |               |
MFADR*LFSTNHKDIGTLYLLFGA*AGVLGTALSLLIRAELGQPGNLLGNDHIYNVIVTAHAFVIIFFIVIPIIIGGFGN*LVPLIIGAPDMAFPRINNISF*LLPPSLLLLLASAIVEAGAGTG*TVYPPLAGNYSHPGASVDLTIFSLHLAGVSSILGAINFITTIINIKPPAITQYQTPLFV*SVLITAVLLLLSLPVLAAGITILLTDRNLNTTFFDPAGGGDPILYQHLF*FFGHPEVYILILPGFGIISHIVTYYSGKKEPFGYIGMV*AMISIGFLGFIV*AHHIFTVGIDVDTRAYFTSATIIIAIPTGVKVFS*LATLHGSNMK*SAAVL*ALGFIFLFTVGGLTGIVLANSSLDIVLHDTYYVVAHFHYVLSIGAVFAIIGGFIH*FPLFSGYTLDQTYAKIHFTIIFIGVNLTFFPQHFLGLSGMPRRYSDYPDAYTT*NILSSVGSFISLTAVILIIFMI*EAFASKRKVLIVEEPSINLE*LYGCPPPYHTFEEPVYIKS


p =

    0.6402

funcs_mt1

score =

  -119


globalAlignment =

MTNIRKSHPLIKIINDAFVDLPAPSNISS**NFGSLLGVCLILQILTGLFLAIHYTPDTLTAFSSVTHICRDVNYG*VIRYIHANGASIFFICLFTHVGRGLYYGSYTFLET*NIGVILLLTTIATAFIGYVLP*GQMSF*GATVITNLLSAIPYIGTDLVE*I*GGFSVDKATLTRFFAFHFILPFVITALAIVHLLFLHETGSNNPTGIPSNADKIPFHPYYTIKDILGILLLITTLLTLTLFAPDLLGDPDNYTPANPLSTPPHIKPE*YFLFAYAILRSIPNKLGGVLALALSILILALIPILHTSKQRSLIFRPLSQCLF*ALIADLLTLT*IGGQPVEHPFIIIGQVASILYFLLILVLMPVAGIIENKLLK*
||:  :||   :|:: :    | |  : :    |: |:: | | : :|| : :|:   ||  :: :|:    : | |  | : :   |  :     |:   :  |    |:    |:||::|:   :|  : :  | : || |   :: | | |   |   : : ||     : :| :  :: : |:: |  : |  | |  :|  : |   :: : |  :   : : |:   |||| || |  |:  || : ::|   :|::    |    :|: | : ::::   :|::: |::  :|  || :  |||:: : |:  :   :  :: |:: |  :    |   :   |:  |  : :   | : ::: ::| : | 
MTH--QSHA-YHIVKPS----P*P--L-T----GA-LSA-L-L-MTSGLAM*FHFHSITLLILGLLTNTL-TI-YQ*-WRDV-TR-EST-YQG-H-HTP-PVQKG----LR---YGIILFITS--EVF--F-FA-G-F-F*-AFYHSS-L-A-P---TP--Q-L-GG-HWPPTGITPLNPLE-V-PLLNT--S-V--L-L-ASGV-SIT--*AHHSLIE-N-NRN-Q-IIQA-LLITILLGL-YFT--LL-QASEYF-ESPFTISDGIYGSTFFV-ATG-FHGLHVIIGSTF-LTIC-FIRQLI-FHFTSKHH-FGFEA-AA*YW-HFV-DVV*LF-L--Y-V-SIY-**GSY-S--FSINSTVNFQLTS-FDN-IQK-


p =

    0.0034

hippo_porf10 = nt2aa(hippo(orfs_valid(4).Start(1) : orfs_valid(4).Stop(1)));
human_porf10 = nt2aa(human(orfs_human(3).Start(4) : orfs_human(3).Stop(4)));

[score, globalAlignment] = nwalign(hippo_porf10, human_porf10)

n = 1000;
globalscores = zeros(n,1);
hippoLen = length(hippo_porf10);
for i = 1:n
    permutedSequence = hippo_porf10(randperm(hippoLen));
    globalscores(i) = nwalign(human_porf10,permutedSequence);
end

figure
buckets = ceil(n/20);
hist(globalscores,buckets)
hold on;
stem(score,1,'c')
title('Determining Alignment Significance: Hippo vs Human Proteins');
xlabel('Score'); ylabel('Number of Sequences');

parmhat = evfit(globalscores);
p = 1 - evcdf(score,parmhat(1),parmhat(2))

score =

   -39


globalAlignment =

DSPVKT*QHSLTQNL*STEQVTPGITAQSYSRVHIDNRVYDLDVGSGHPNGAAAIKGSFVQRLKSYVI*VQTGAIQVSFYLLYISPSTKGQEK*GLLPISALKLINDIVLT*LNSINIPALDQGTVAMAEPGNCIKLKPLHQRFKSSSQQNVYHQHPYTCRTHPSSHSIPNTSRTKNP
 | :      |: :| : | :  |:   :   |:  |:: |:| : :: |    :   ::: : :::|    | | |  |  ::  |    |  |:: :| | :|  :: | ||: :: | :   :|:   : |     |   : | |  |:|    |   :    |:: |    :  
MSSLLL-IIILALSL-AYE*LQKGLD-*A-ELVYSLNKTNDFD-SLNYDNHIYQMPLIYINIILAFTI-SLLG-ILV--YRSHLMSSLLCLE--GII-LS-LFIIATLI-T-LNTHSLLA-NIVPIAILVFAACEAAVGL-ALLVSIS--NTYGL-DY---V----HNL-NL--LQC-


p =

    0.0132
    
 %% funcs mt1.m

hippo_porf10 = nt2aa(hippo(orfs_valid(4).Start(1) : orfs_valid(4).Stop(1)));
human_porf10 = nt2aa(human(orfs_human(3).Start(4) : orfs_human(3).Stop(4)));

[score, globalAlignment] = nwalign(hippo_porf10, human_porf10)

n = 1000;
globalscores = zeros(n,1);
hippoLen = length(hippo_porf10);
for i = 1:n
    permutedSequence = hippo_porf10(randperm(hippoLen));
    globalscores(i) = nwalign(human_porf10,permutedSequence);
end

figure
buckets = ceil(n/20);
hist(globalscores,buckets)
hold on;
stem(score,1,'c')
title('Determining Alignment Significance: Hippo vs Human Proteins');
xlabel('Score'); ylabel('Number of Sequences');

parmhat = evfit(globalscores);
p = 1 - evcdf(score,parmhat(1),parmhat(2))

%% end of script

[score, localAlignment] = swalign(hippoProtein, humanProtein)

n = 50;
localscores = zeros(n,1);
humanLen = length(humanProtein);
aminos = aacount(humanProtein);
for i = 1:n
    randSequence = randseq(humanLen, 'FROMSTRUCTURE', aminos);
    localscores(i) = swalign(hippoProtein,randSequence);
end

figure
buckets = ceil(n/5);
hist(localscores,buckets)
hold on;
stem(score,3,'c')
title('Monte Carlo Simulation');
xlabel('Score'); ylabel('Number of Sequences');
parmhat = evfit(localscores);
x = min(localscores):max([localscores; score]);
y = evpdf(x,parmhat(1),parmhat(2));
[v, c] = hist(localscores,buckets);
binWidth = c(2) - c(1);
scaleFactor = n*binWidth;
plot(x,scaleFactor*y,'r');
hold off;


p = 1 - evcdf(score,parmhat(1),parmhat(2))

score =

   1.4023e+03


localAlignment =

*HATTRHINMIYH-HPIHISDPIYYLSTENLKTHLPPKP*DYSSHNTKTAYPLRNEMNENLFASFITPTILGLPLVTLIIMFPSILFPAPTRLITNRLVSIQQ*LIQLVSKQIMNIHNHKGQT*TLILISLILFIGSTNLLGLLPHSFTPTTQLSINLGIAIPL*AGTVIIGFRNKTKISLAHFLPQGTPTPLIPMLVIIETISLFIQPIALAVRLTANITAGHLLMHLIGGATLALINISITTALITFIILVLLTALEFAVAIIQAYVFTLLVSLYLHDNT**PTKPTHTI**TQVPDLLQEPSQPY**RRA*PYDSTLTP--LSY*RQD*LPIS*QYISDDEM*SEKAPFKATTHQSYKKDFATE*SYLLSPKSYF-SQASSEPFTTQASLLLLN*ADVDHPQASTL*TH*RCHF*TPPFY*PLVSPLPEPTTV**KAIENKYSKPSSSQSP*VCTSHYCKLQNTMKPPLQSQMGFTAQLSL*PQAFMDYM*LLAPLS*LYASYAN*NSTSRQITTLASRPPPDTDTS*M*SDYSSTCPFIDEVHSSFSIKTVQLTSNQLASVHSGKEQLT***HY*QTPH*PLYWSSLPSDSHN*TPTQKKQAPM-NADLTL*DQPAYLSL*NSS**PSHSFFST*RSPSYFLSHGQP-KQQT*KPYSL*PSP*SHS*QSA*PTNELKRD*NEPN-MVFSLKQNK*FRLIKL*TNS*IPSVFSIYKYYHGLYNI-PR-RTVNISI-PLNILTPMSRRNDIITIY-HSN-SHH-P-KCTL-HPSQHNANYSTSFRSM*S-SPRTIATSNGIKHIRYRLRTKPKPSPMLKYIIPTIILIPLT*ISKNSII*TNTTAHSLLISFTSLLLLNQFNDNSLNFSPMFFSDPLSTPLLILTIWLLPLILIASQ-SHLLK-EPPTRKKL-FITILVT-LQTFLIITFSAIELILFYILFEATLIPTL-I-IITR*GNQT-ERLNAGLYFLFYTLIGSLPLLVALIYIQNITGSLNFLILQY*TQAVSNS*SNVFL*LACI-IAFMVKIPLYGLHL*LPKAHVEAPIAGSIVLAAILLKLGGYGILRITTILNPLTEIMAYPFIILSL*GIIMTSSICLRQTDLKSLIAYSSVSHIALVIVAILIQTP*SYIGATALIIAHGLTSS-ILFCLANSNYE-RIHSRT--IILARGLQTLLPLIAA**LLASLTNLALPPSINLVGELLVIMSSFS*SNITIILIGT-NIIITALYSLYILTTTQRGKYTHHINNITPSFTRENALI-A-L-HILPLLLLSLNPKIILGPLYCKHSLRKTLDCESTNRSSNPSYLPRKHARTANSCPHI*QYGFPRLLKDGSYPLVLGTKKLVQLQIKAINLFSS-TTLTILFVLTLPIIITNTNIYKSDKYPTYVKNTVSSAFLISLVPIIAFTNTGQEIIISN*H*ITIQTLKLTLSFKADYFSIVFAPVALFVTWSIMEFSIWYIHSDPHINQFFKYLLLFLITIIVLVTANNLFQLFIG*EGVGIMSFLLIG**HGRTDANTAAIQAILYNRIGDVGFIIAIA*FLSNLNT*DIQQIFIINPTHSNLPLIGLILAATGKSAQFGLHP*LPSAMEGPTPVSALLHSSTIVVAGVFLLIRFYPMIENNNLIQTITICLGAITTLFTAICALTQNDIKKIIAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFICSGSIIHNLNNEQDIRKIGGLFKTIPFTTTTLIVGSIALTGVPFLTGFYSKDLIIEAANTSYSNA*ALLITLVATSLTAVYSTRIIFFALLGHPRFPTSTLINENNPLLLNSLKRLIAGSIFAGFILSHNLPPITTPL-ITIPPYLKMTALAVTILGFTLAFEITLNTQNLKHKHPTNSFKFSTLLGYFPTIIHRLPPHLSLTASQKLASSLLDSA*LENILPKSIAHAQLKLSTLVSNQKGLIKIYFLSFLITIPLSLILFN-P
|:| |::  | :| :| |    | :  |:|:| :|||     : :| |    ||::||||||||||:||||||| ::|||:|| :|:|:   ||:|||:: |||||:|:|||:::||| ||:||:|||:|||:||::|||||||||||||||||||||::|||||||:||||||:| | :|||||||||||||||:|||||||||:|||||||||||||||||||||||||::|||: :|:: ::|| | ||:||| ||:|||:||||||||||||||||||||||:    |||: : |  | |||| ||  ||| | | ||   |||   | |  | | :|  || |||  |||||:  || | | ||||| |:::| ||  |||||| |  |  ||  :   |||  |  |: | || |:|       || | ||||: |:|  |  | |   |      | |:| : |     : |||  |||||  |:  ||| || | || |||  | :   ||||:|||||   ||  ||    | |||   :: :  | |||| |: :   ||  |  |  |:|| || |   |  | || | |  |: |: :|  ||   || || ||||||   |::  |: ||::  | :| || | :     || || ||  : | :: |:||:||| ::   |:    : | :   :  || |:   | |   ::      :: |  |         | ::  :   |: :   |  | | : |   |:   :   :|| ||  :| :: |    |   :| :  |  :||  | : :  |   :| :   |    :     ::::  || :   | | | :|  | ::|    | :| |     :: : |::    |:| | : | | | | :: | :  : : || |:| : |  | |      |   :| || :  | :|  : : | |   : :: |:   : :     :  ::::   |:  | :: :  :|: |       |  | :  ::|           |:          |  :|  :  : | :  ::      ::: |  |: :  ||:  : :      :   : | : :       :|:   : ||: :  |: : | ||::    |:|  ::  |  :    :::     | |    :|: |   :     :       | :   : :  |:|  :|:   :  | :: |   ||::  : : : |  |:     : :::|  :  ::::|  |   : :| ::| ||||| :||||||||:||:||  : :  :     : : |||||||:|:: :: |||| |  |  ||: | :|  |:::|| |||: |:|:|:||| |   |    ||:||||||| | || :|:|||| |||||:| ||||||||||:|||: ||:|||:|||||||||:|||||::|||||||||||||||||||:|||||:||::|:|||||||||||||||||:|||:|:|||: : |:|| ||| ::| : |  ||:||:|||:|||||:|||||||||:|||||||||||||||||||:||||||:|: ||: ||||:|:|||||||||:|:|||||||||||:||||||||||||||||||||||||||||||||||||||:|||||||||||||||||||||:||||:|:|:| :||:||:|:|||||||||| |||:|| ||:||||| |||:|||||::||||||:::| |:||||| | |||||| ||| :||| |||:||||::::|: | ::|:  ||| |||:||||||:||:  |::::  |::|| | |  :| ||::||::|:| ||  |:|:| :||:|   ||| :|||::|||:|:: |:: | ::|:|||:||:|||||:: : |:|:|:: |
*NAPTKYYRMAHHNYP-HTPYTIPHHPTKNIKHKLPPTSLTKAHKNKKL*QTLRTKMNENLFASFIAPTILGLPAAVLIILFPPLLIPTSKYLINNRLITTQQ*LIKLTSKQMITIHNTKGRT*SLILVSLIIFIATTNLLGLLPHSFTPTTQLSINLAMAIPL*AGAVIIGFRSKIKNALAHFLPQGTPTPLIPILVIIETISLLIQPIALAVRLTANITAGHLLMHLIGSTTLAISTINLPSTLIIFTILILLTILEIAVALIQAYVFTLLVSLYLHDNT**PTNHMPII**NPAHDP*QGPSQPS**PPA*PCDFTSTP*RSSY*AY-*-PTH*PYTNDGAM*HEKAHTKATTHHLSKKAFDTG*SYLL-PQKFFSSQDFSEPFTTPA*PLPPN*EGTGPQQASPR*IP*KSHS*THPYYSHQEYQSPELTIV**KTTETK*FKHCSLQFYWVSILPSYKPQSTSSLPSPFPTASTAQHFL*PQASTDFTSLLAQLSSLSASSAN*YFTLHPNITLASKPPPDTGIL*MWFDYFCMSPSIDE-GLTLLV*IVPLTSN*LVLTTFKKE**TSP*F**STPS*PYY**LLHFDYHNSTAT-*KNPPLTSAASTLYPPPASLSP*NSS**LLPSYYLI*KLPSFYPYH-EPYKQLTCH**LCHPSY*SSS*P*VWPMSDYKKD*TEPNWYIV*TKRMISTH*IMIIIFTKCPS-FT*ILY*HLPSHF*EY*YIAHTSCPPYYA*KE*YYRCSL*LLS*PSTPTPS*PILCLLPY*SLPPAKQRWA-*PY*SQSPTHMA*TTYIT-*TYS-NAKTN-RPN-NYITTTDMTFQKT---HN-LNQHNHPQPNY*HHPSTIF*PNQQQPI*L-F-PNLFLRPPNNP-PPNTNYLTPTPHNHGKPTPLIQ*TTITKKTLPLYTNLPTNLLNYNIHSHRT-NHIL-YLLRNHT-YPHLGYHHPMRQPARTPERRHI-LPILHPSRLPS-PTHRTNLHSQHPRLTKHSTTHPHCPRTIKLL-SQQ-LNMTSLHNSFYSKDTSLRTPLMTP*SPCRSPHRWVNSTCRSTLETRRLWHNTPHTHSQPPDKTHSLP-LPCTIPMRHNYNKLHLPTTN-RPKIAHCILFNQPHSPRSNSHSHP-NPLKLHRRSHSHNRPRTYILITILPSKLKLRTHSQSHHNPLSRTSNS-TPTNSFLMTFSKPR*PRLTPHY*PTGRTLCASNHVLLIKYHSPTYRTQHTSHSPILPLHIYHNTMGLTHPPH*QHKTLIHTRKHPHVHTPIPHSPPIPQPRHHYRVFLL*I*FNQNIR--L*I*Q-QRLTTP-YLPRKLTRTANSCPHV*QHGFLNF*RITAIHWS*APRILVQLQIKVITMHTTITTLT-LTSLIPPILTTLVNPNKKNSYPHYVKSIVASTFIISLFPTTIFMCLDQEVIISN*H*ATTQTTQLSLSFKLDYFSIIFIPVALFVTWSIIEFSL*YINSDPNINQFFKYLLIFLITILILVTANNLFQLFIG*EGVGIISFLLIS**YARADANTAAIQAILYNRIGDIGFILALA*FILHSNS*DPQQIALLNANPSLTPLLGLLLAAAGKSAQLGLHP*LPSAIEGPTPVSALLHSSTIVVAGIFLLIRFHPLAENSPLIQTLTLCLGAITTLFAAVCALTQNDIKKIVAFSTSSQLGLIIVTIGINQPHLAFLHICTHAFFKAILFMCSGSIIHNLNNEQDIRKIGGLLKTIPLTSTSLTIGSLALAGIPFLTGFYSKDHIIETANISYTNA*ALSITLIATSLTSAYSTRIILLTLTGQPRFPTLTNINENNPTLLNPIKRLAAGSLFAGFLITNNISP-ASPFQTTIPLYLKLTALAVTFLGLLTALDLNYLTNKLKIKSPLCTFYFSNILGFYPSITHRTIPYLGLLTSQNLPLLLLDLT*LEKLLPKTISQHQISTSIITSTQKGIIKLYFLSFFFPLILTLLLIT*P


p =

     0

 

