% Retrieve sequence from genbank
WM = getgenbank('AP008987')

% Report basic statistical descriptions of retrieved sequence
mammoth_full = fastaread('~/Downloads/mammoth.fasta')
mammoth = mammoth_full.Sequence
length(mammoth)
basecount(mammoth)
ntdensity(mammoth)

figure
dimers = dimercount(mammoth, 'chart', 'bar')
title('Dimer Count for mammuthus primigenius DNA sequence')


mammothORFs = seqshoworfs(mammoth, 'geneticcode', 2, 'frames', [1,2,3,-1,-2,-3])


% Find protein-coding genes and translate two of them
features = featureparse(WM, 'Sequence', true)
coding_seqs = features.CDS;
coding_seqs_id = sprintf('%s ', coding_seqs.gene)

WM_CYTB = coding_seqs(13).Sequence
WM_COX1 = coding_seqs(3).Sequence

% Visualise complete genome with protein coding regions
[h, l] = featureview(WM, {'CDS', 'tRNA', 'rRNA', 'D_loop'}, [2 1 2 2 2], 'Fontsize', 9);
legend(h, l, 'interpreter', 'none');
title('Mammuthus primigenius mitochondria DNA, complete genome')

% Translate protein coding regions
WM_CYTBaa = nt2aa(WM_CYTB, 'GeneticCode', 'Vertebrate Mitochondrial');
disp(seqdisp(WM_CYTBaa))

WM_COX1aa = nt2aa(WM_COX1, 'GeneticCode', 'Vertebrate Mitochondrial');
disp(seqdisp(WM_COX1aa))

% Visualise amino acid counts
figure
subplot(2,1,1)
CYTBaaCount = aacount(WM_CYTBaa, 'chart', 'bar');
title('Histogram of Amino Acid Count for CYTB protein');

figure
subplot(2,1,2)
COX1aaCount = aacount(WM_COX1aa, 'chart', 'bar');
title('Histogram of Amino Acid Count for COX1 protein');

%% Plot phylogenetic trees for first 5 species vs mammoth
data = {'Columbian_Mammoth'      'NC_015529';
        'African_Sav_Elephant'     'NC_000934';
        'African_Forest_Elephant'          'JN673263'  ;
        'Asiatic_Elephant' 'EF588275';
        'American_Mastodon'       'EF632344'; 
        'Woolly_Mammoth'    'AP008987';
       };
        
for ind = 1:6
    seqs(ind).Header   = data{ind,1};
    seqs(ind).Sequence = getgenbank(data{ind,2},...
                                    'sequenceonly', true);
end

% UPGMA tree
distances = seqpdist(seqs,'Method','Jukes-Cantor','Alphabet','NT');
tree = seqlinkage(distances,'UPGMA',seqs)
h = plot(tree,'orient','top');
ylabel('Genetic Distance')
title('UPGMA Phylogenetic Tree')
set(h.terminalNodeLabels,'Rotation',65)
        
% Neighbour-joining tree
distances = seqpdist(seqs,'Method','Jukes-Cantor','Alphabet','DNA')
njtree = seqneighjoin(distances,'equivar',seqs)
h = plot(tree,'orient','top');
ylabel('Genetic Distance')
title('Neighbour Joining Phylogenetic Tree')
set(h.terminalNodeLabels,'Rotation',65)


%% Get diversified species genomes and extract CYTB and COX1 proteins

% Columbian Mammoth
CM = getgenbank('NC_015529')
features = featureparse(CM, 'Sequence', true);
coding_seqs = features.CDS;
CM_CYTB = coding_seqs(13).Sequence;
CM_COX1 = coding_seqs(3).Sequence;
CM_CYTBaa = nt2aa(CM_CYTB, 'GeneticCode', 'Vertebrate Mitochondrial');
CM_COX1aa = nt2aa(CM_COX1, 'GeneticCode', 'Vertebrate Mitochondrial');

% African Savanna Elephant
ASE = getgenbank('NC000934')
features = featureparse(ASE, 'Sequence', true);
coding_seqs = features.CDS;
ASE_CYTB = coding_seqs(13).Sequence;
ASE_COX1 = coding_seqs(3).Sequence;
ASE_CYTBaa = nt2aa(ASE_CYTB, 'GeneticCode', 'Vertebrate Mitochondrial');
ASE_COX1aa = nt2aa(ASE_COX1, 'GeneticCode', 'Vertebrate Mitochondrial');

% Asiatic Elephant
AE = getgenbank('EF588275')
features = featureparse(AE, 'Sequence', true);
coding_seqs = features.CDS;
AE_CYTB = coding_seqs(13).Sequence;
AE_COX1 = coding_seqs(3).Sequence;
AE_CYTBaa = nt2aa(AE_CYTB, 'GeneticCode', 'Vertebrate Mitochondrial');
AE_COX1aa = nt2aa(AE_COX1, 'GeneticCode', 'Vertebrate Mitochondrial');

% American Mastodon
AM = getgenbank('EF632344')
features = featureparse(AM, 'Sequence', true);
coding_seqs = features.CDS;
AM_CYTB = coding_seqs(13).Sequence;
AM_COX1 = coding_seqs(3).Sequence;
AM_CYTBaa = nt2aa(AM_CYTB, 'GeneticCode', 'Vertebrate Mitochondrial');
AM_COX1aa = nt2aa(AM_COX1, 'GeneticCode', 'Vertebrate Mitochondrial');

% African White Lion
AWL = getgenbank('KF907306')
features = featureparse(AWL, 'Sequence', true);
coding_seqs = features.CDS;
AWL_CYTB = coding_seqs(13).Sequence;
AWL_COX1 = coding_seqs(3).Sequence;
AWL_CYTBaa = nt2aa(AWL_CYTB, 'GeneticCode', 'Vertebrate Mitochondrial');
AWL_COX1aa = nt2aa(AWL_COX1, 'GeneticCode', 'Vertebrate Mitochondrial');


%% Compute genetic distance between each pair seperately, for each type of protein coding sequence

% Create structures for NT sequences
data_COX1 = {'Columbian_Mammoth'      CM_COX1;
        'African_Sav_Elephant'     ASE_COX1;
        'African_White_Lion'       AWL_COX1;
        'Asiatic_Elephant' AE_COX1;
        'American_Mastodon'       AM_COX1; 
        'Woolly_Mammoth'    WM_COX1;
       };
for ind = 1:6
    seqs_COX1(ind).Header   = data_COX1{ind,1};
    seqs_COX1(ind).Sequence = data_COX1{ind,2};
end


data_CYTB = {'Columbian_Mammoth'      CM_CYTB;
        'African_Sav_Elephant'     ASE_CYTB;
        'African_White_Lion'       AWL_CYTB;
        'Asiatic_Elephant' AE_CYTB;
        'American_Mastodon'       AM_CYTB; 
        'Woolly_Mammoth'    WM_CYTB;
       };
for ind = 1:6
    seqs_CYTB(ind).Header   = data_CYTB{ind,1};
    seqs_CYTB(ind).Sequence = data_CYTB{ind,2};
end


% Create structures for AA sequences
data_CYTBaa = {'Columbian_Mammoth'      CM_CYTBaa;
        'African_Sav_Elephant'     ASE_CYTBaa;
        'African_White_Lion'       AWL_CYTBaa;
        'Asiatic_Elephant' AE_CYTBaa;
        'American_Mastodon'       AM_CYTBaa; 
        'Woolly_Mammoth'    WM_CYTBaa;
       };
for ind = 1:6
    seqs_CYTBaa(ind).Header   = data_CYTBaa{ind,1};
    seqs_CYTBaa(ind).Sequence = data_CYTBaa{ind,2};
end

data_COX1aa = {'Columbian_Mammoth'      CM_COX1aa;
        'African_Sav_Elephant'     ASE_COX1aa;
        'African_White_Lion'       AWL_COX1aa;
        'Asiatic_Elephant' AE_COX1aa;
        'American_Mastodon'       AM_COX1aa; 
        'Woolly_Mammoth'    WM_COX1aa;
       };
for ind = 1:6
    seqs_COX1aa(ind).Header   = data_COX1aa{ind,1};
    seqs_COX1aa(ind).Sequence = data_COX1aa{ind,2};
end

% Compute genetic distance between each pair for NT sequences
distances_CYTB = seqpdist(seqs_CYTB, 'Method', 'Jukes-Cantor', 'Alphabet', 'NT', 'PairwiseAlignment', true);
D = squareform(distances_CYTB);
imagesc(D);
legend
colorbar
title('Pairwise distances - CYTB NT') 

distances_COX1 = seqpdist(seqs_COX1, 'Method', 'Jukes-Cantor', 'Alphabet', 'NT', 'PairwiseAlignment', true);
D = squareform(distances_COX1);
imagesc(D);
legend
colorbar
title('Pairwise distances - COX1 NT') 

% Compute genetic distance between each pair for AA sequences
distances_CYTBaa = seqpdist(seqs_CYTBaa, 'Method', 'Jukes-Cantor', 'Alphabet', 'AA', 'PairwiseAlignment', true);
D = squareform(distances_COX1);
imagesc(D);
legend
colorbar
title('Pairwise distances - CYTB AA')

distances_COX1aa = seqpdist(seqs_COX1aa, 'Method', 'Jukes-Cantor', 'Alphabet', 'AA', 'PairwiseAlignment', true);
imagesc(D);
legend
colorbar
title('Pairwise distances - COX1 AA')


%% Generate rooted phylogenetic trees using an outgroup

% Prepare outgroup data
GGC = getgenbank('AY235571')
features = featureparse(GGC, 'Sequence', true);
coding_seqs = features.CDS;
coding_seqs_id = sprintf('%s ', coding_seqs.gene)
GGC_CYTB = coding_seqs(11).Sequence
GGC_COX1 = coding_seqs(3).Sequence
GGC_CYTBaa = nt2aa(GGC_CYTB, 'GeneticCode', 'Vertebrate Mitochondrial')
GGC_COX1aa = nt2aa(GGC_COX1, 'GeneticCode', 'Vertebrate Mitochondrial')

% Add outgroup to data and plot phylogenetic trees
% Tree 1
data_CYTB = {'Columbian_Mammoth'      CM_CYTB;
        'African_Sav_Elephant'     ASE_CYTB;
        'African_White_Lion'       AWL_CYTB;
        'Asiatic_Elephant' AE_CYTB;
        'American_Mastodon'       AM_CYTB; 
        'Woolly_Mammoth'    WM_CYTB;
        'Chicken'       GGC_CYTB;
       };
for ind = 1:7
    seqs_CYTB(ind).Header   = data_CYTB{ind,1};
    seqs_CYTB(ind).Sequence = data_CYTB{ind,2};
end
distances = seqpdist(seqs_CYTB,'Method','Jukes-Cantor','Alphabet','NT', 'PairwiseAlignment', true');
tree = seqlinkage(distances,'UPGMA',seqs_CYTB);
h = plot(tree,'orient','top');
ylabel('Genetic Distance')
title('Phylogenetic Tree - CYTB NT')
set(h.terminalNodeLabels,'Rotation',65)

% Tree2
data_COX1 = {'Columbian_Mammoth'      CM_COX1;
        'African_Sav_Elephant'     ASE_COX1;
        'African_White_Lion'       AWL_COX1;
        'Asiatic_Elephant' AE_COX1;
        'American_Mastodon'       AM_COX1; 
        'Woolly_Mammoth'    WM_COX1;
        'Chicken'       GGC_COX1;
       };
for ind = 1:7
    seqs_COX1(ind).Header   = data_COX1{ind,1};
    seqs_COX1(ind).Sequence = data_COX1{ind,2};
end
distances = seqpdist(seqs_COX1,'Method','Jukes-Cantor','Alphabet','NT', 'PairwiseAlignment', true');
tree = seqlinkage(distances,'UPGMA',seqs_COX1);
h = plot(tree,'orient','top');
ylabel('Genetic Distance')
title('Phylogenetic Tree - COX1 NT')
set(h.terminalNodeLabels,'Rotation',65)

% Tree 3
data_COX1aa = {'Columbian_Mammoth'      CM_COX1aa;
        'African_Sav_Elephant'     ASE_COX1aa;
        'African_White_Lion'       AWL_COX1aa;
        'Asiatic_Elephant' AE_COX1aa;
        'American_Mastodon'       AM_COX1aa; 
        'Woolly_Mammoth'    WM_COX1aa;
        'Chicken'       GGC_COX1aa;
       };
for ind = 1:7
    seqs_COX1aa(ind).Header   = data_COX1aa{ind,1};
    seqs_COX1aa(ind).Sequence = data_COX1aa{ind,2};
end
distances = seqpdist(seqs_COX1aa,'Method','Jukes-Cantor','Alphabet','AA', 'PairwiseAlignment', true');
tree = seqlinkage(distances,'UPGMA',seqs_COX1aa);
h = plot(tree,'orient','top');
ylabel('Genetic Distance')
title('Phylogenetic Tree - COX1 AA')
set(h.terminalNodeLabels,'Rotation',65)

% Tree 4
data_CYTBaa = {'Columbian_Mammoth'      CM_CYTBaa;
        'African_Sav_Elephant'     ASE_CYTBaa;
        'African_White_Lion'       AWL_CYTBaa;
        'Asiatic_Elephant' AE_CYTBaa;
        'American_Mastodon'       AM_CYTBaa; 
        'Woolly_Mammoth'    WM_CYTBaa;
        'Chicken'       GGC_CYTBaa;
       };
for ind = 1:7
    seqs_CYTBaa(ind).Header   = data_CYTBaa{ind,1};
    seqs_CYTBaa(ind).Sequence = data_CYTBaa{ind,2};
end
distances = seqpdist(seqs_CYTBaa,'Method','Jukes-Cantor','Alphabet','AA', 'PairwiseAlignment', true');
tree = seqlinkage(distances,'UPGMA',seqs_CYTBaa);
h = plot(tree,'orient','top');
ylabel('Genetic Distance')
title('Phylogenetic Tree - CYTB AA')
set(h.terminalNodeLabels,'Rotation',65)


%% Perform multiple alignments

% CYTB NT alignment
ma = multialign(seqs_CYTB)
seqalignviewer(ma)

% COX1 NT alignment
ma2 = multialign(seqs_COX1)
seqalignviewer(ma2)




