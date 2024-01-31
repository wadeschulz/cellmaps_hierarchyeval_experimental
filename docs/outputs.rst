=======
Outputs
=======

The **cellmaps_hierarchyeval** tool produces several output files within the specified directory. These files are
the results of the enrichment analysis, as well as, logs and essential metadata.

Enriched Hierarchy and Parent Network
-------------------------------------
The enriched hierarchy and parent network are provided in CX2 format:

- ``hierarchy.cx2``:
    This is the enriched hierarchy network file that integrates the results of the enrichment analysis into the hierarchy, formatted in CX2.

Example of annotated node in the hierarchy ("HCX::members" value was removed for readability):

.. code-block:: json

    {
        "id": 226,
        "v": {
            "CD_MemberList": "LRRFIP2 CNN3 SEPTIN5 TNNC1 SEPTIN7 FAM216A GPX8 PRKRIP1 ACTN4 SPRYD3 LSM6 CDC42EP4 SPECC1L BZW2 FRMD1 HTRA1 SZT2 BBOX1 BRICD5 MYH9 PDRG1 TPM3 RAI14 LIMCH1 CTPS1 SIPA1L1 SEPTIN9 NEXN APPL1 LUZP1 WASHC3 PPP1R12A SEPTIN3 SEPTIN10 GABRA3 TAX1BP3 NCOA5 GSN MAP2 ATP6V1H DMWD",
            "CD_MemberList_Size": 41,
            "CD_MemberList_LogSize": 5.358,
            "CD_AnnotatedMembers": "",
            "CD_AnnotatedMembers_Size": 0,
            "CD_AnnotatedMembers_Overlap": 0.0,
            "CD_AnnotatedMembers_Pvalue": 0.0,
            "HiDeF_persistence": 81,
            "CD_CommunityName": "C5044",
            "CD_Labeled": true,
            "HCX::isRoot": false,
            "HCX::members": [],
            "CORUM_terms": "",
            "CORUM_FDRs": "",
            "CORUM_jaccard_indexes": "",
            "CORUM_overlap_genes": "",
            "GO_CC_terms": "GO:0031105|GO:0005940|GO:0032156",
            "GO_CC_descriptions": "septin complex|septin ring|septin cytoskeleton",
            "GO_CC_FDRs": "3.150973655449002e-07|3.150973655449002e-07|6.709368907329383e-07",
            "GO_CC_jaccard_indexes": "0.11904761904761904|0.11904761904761904|0.11627906976744186",
            "GO_CC_overlap_genes": "SEPTIN5,SEPTIN3,SEPTIN10,SEPTIN9,SEPTIN7|SEPTIN5,SEPTIN3,SEPTIN10,SEPTIN9,SEPTIN7|SEPTIN5,SEPTIN3,SEPTIN10,SEPTIN9,SEPTIN7",
            "HPA_terms": "Actin filaments",
            "HPA_FDRs": "6.450286324867403e-58",
            "HPA_jaccard_indexes": "0.375",
            "HPA_overlap_genes": "TPM3,BRICD5,CTPS1,SEPTIN5,SEPTIN3,GPX8,MYH9,SEPTIN10,CNN3,LUZP1,BBOX1,SPECC1L,PRKRIP1,LSM6,SEPTIN7,PPP1R12A,BZW2,LRRFIP2,LIMCH1,FRMD1,CDC42EP4,DMWD,NCOA5,PDRG1,FAM216A,SIPA1L1,NEXN,SZT2,TNNC1,SPRYD3,ATP6V1H,SEPTIN9,GSN,RAI14,ACTN4,TAX1BP3",
            "n": "C5044",
            "r": "C5044"
        }
    }

- ``hierarchy_parent.cx2``:
    The reference parent network from which the hierarchy was generated, formatted in CX2. Copy from input.

Node Attributes
---------------
The attributes for each node within the enriched hierarchy:

- ``hierarchy_node_attributes.tsv``:
    A TSV file containing attributes for each node, which includes information such as enriched terms, their descriptions, and related statistical data.

.. code-block::

    name	represents	CD_MemberList	CD_MemberList_Size	CD_MemberList_LogSize	CD_AnnotatedMembers	CD_AnnotatedMembers_Size	CD_AnnotatedMembers_Overlap	CD_AnnotatedMembers_Pvalue	HiDeF_persistence	CD_CommunityName	CD_Labeled	HCX::isRoot	HCX::members	CORUM_terms	CORUM_FDRs	CORUM_jaccard_indexes	CORUM_overlap_genes	GO_CC_terms	GO_CC_descriptions	GO_CC_FDRs	GO_CC_jaccard_indexes	GO_CC_overlap_genes	HPA_terms	HPA_FDRs	HPA_jaccard_indexes	HPA_overlap_genes
    C5044	C5044	LRRFIP2 CNN3 SEPTIN5 TNNC1 SEPTIN7 FAM216A GPX8 PRKRIP1 ACTN4 SPRYD3 LSM6 CDC42EP4 SPECC1L BZW2 FRMD1 HTRA1 SZT2 BBOX1 BRICD5 MYH9 PDRG1 TPM3 RAI14 LIMCH1 CTPS1 SIPA1L1 SEPTIN9 NEXN APPL1 LUZP1 WASHC3 PPP1R12A SEPTIN3 SEPTIN10 GABRA3 TAX1BP3 NCOA5 GSN MAP2 ATP6V1H DMWD	41	5.358		0	0	0	81	C5044	TRUE	FALSE	[4002, 92, 4446, 3572, 36, 2324, 4131, 3546, 1008, 294, 3722, 4786, 1923, 4241, 4756, 2307, 4804, 4970, 2326, 35, 1009, 4110, 633, 4169, 2733, 4858, 4775, 4963, 2368, 287, 1215, 4440, 3016, 2986, 4927, 290, 3566, 632, 1033, 289, 4262]					GO:0031105|GO:0005940|GO:0032156	septin complex|septin ring|septin cytoskeleton	3.150973655449002e-07|3.150973655449002e-07|6.709368907329383e-07	0.11904761904761904|0.11904761904761904|0.11627906976744186	SEPTIN5,SEPTIN3,SEPTIN10,SEPTIN9,SEPTIN7|SEPTIN5,SEPTIN3,SEPTIN10,SEPTIN9,SEPTIN7|SEPTIN5,SEPTIN3,SEPTIN10,SEPTIN9,SEPTIN7	Actin filaments	6.45E-58	0.375	TPM3,BRICD5,CTPS1,SEPTIN5,SEPTIN3,GPX8,MYH9,SEPTIN10,CNN3,LUZP1,BBOX1,SPECC1L,PRKRIP1,LSM6,SEPTIN7,PPP1R12A,BZW2,LRRFIP2,LIMCH1,FRMD1,CDC42EP4,DMWD,NCOA5,PDRG1,FAM216A,SIPA1L1,NEXN,SZT2,TNNC1,SPRYD3,ATP6V1H,SEPTIN9,GSN,RAI14,ACTN4,TAX1BP3
    C5285	C5285	NAA15 NAA16 NAA50 HYPK	4	2		0	0	0	19	C5285	TRUE	FALSE	[2258, 2257, 2565, 4598]					GO:0031415|GO:0031414	NatA complex|N-terminal protein acetyltransferase complex	1.943122029855394e-07|2.2862713150320575e-06	0.6|0.3333333333333333	NAA15,NAA16,NAA50|NAA15,NAA16,NAA50

Logs and Metadata
-----------------
- ``error.log``:
    Contains error messages and exceptions that might have occurred during execution.

- ``output.log``:
    Provides detailed logs about the steps performed and their outcomes.

- ``ro-crate-metadata.json``:
    Metadata in RO-Crate_ format, a community effort to establish a lightweight approach to packaging research data with their metadata.

    It contains general information about the data i.a. ID, Type, Name, Description, contextual definitions,
    Software detail, as well as datasets details of each individual part of the data.

    For example, the metadata for the content of hierarchy.cx provides unique id, context, type, url, name, keywords, etc.
    The url can be used to view the hierarchy in Cytoscape_ Web.

    .. code-block:: json

        {
          "@id": "00000000-0000-0000-0000-000000000000:dataset::4.hierarchy",
          "@context": {
            "@vocab": "https://schema.org/",
            "evi": "https://w3id.org/EVI#"
          },
          "metadataType": "https://w3id.org/EVI#Dataset",
          "url": "https://idekerlab.ndexbio.org/cytoscape/network/aaaaaaaa-aaaa-aaaa-aaaa-aaaaaaaaaaaa",
          "name": "Output Dataset",
          "keywords": [
            "CM4AI",
            "Example",
            "interactome",
            "ppi",
            "network",
            "CX2",
            "file",
            "hierarchy",
            "network",
            "HCX"
          ],
          "description": "CM4AI Example Example input dataset AP-MS edgelist download|IF microscopy merged embedding AP-MS edgelist download|IF microscopy Example input dataset hierarchy model Hierarchy network file",
          "author": "cellmaps_generate_hierarchy",
          "datePublished": "2023-09-21",
          "version": "0.1.0a11",
          "associatedPublication": null,
          "additionalDocumentation": null,
          "format": "HCX",
          "schema": {},
          "generatedBy": [],
          "derivedFrom": [],
          "usedBy": [],
          "contentUrl": "path/hierarchy.hcx"
        }

    Additionally, it contains Computation Details, name, description, Run By etc.

.. _RO-Crate: https://www.researchobject.org/ro-crate
.. _Cytoscape: https://cytoscape.org/
