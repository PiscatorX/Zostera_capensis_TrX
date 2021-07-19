# *Zostera capensis* transcriptome assembly

## Ion Torrent Sequencing

* Library sequenced on [Ion S5](https://www.thermofisher.com/order/catalog/product/A27212#/A27212) with [Ion 540 Kit-Chef](https://www.thermofisher.com/order/catalog/product/A30011?SID=srch-srp-A30011#/A30011?SID=srch-srp-A30011)
* Reads generate using [IonXpressRNA](https://www.thermofisher.com/order/catalog/product/4475485#/4475485) Barcodes.
* Background reading [blogpost](https://labmedicineblog.com/2018/01/05/template-preparation-clonal-amplification-and-cluster-generation-oh-my-step-two-in-an-ngs-setup/)

## Raw read data

* Sequencing reads were received as unmapped BAM (uBAM) files, reported [here](http://147.47.77.94/ion-docs/Technical-Note---Transition-from-SFF-to-BAM-format_37421247.html).
There is biostars [thread](https://www.biostars.org/p/138937/) on it
* Read converted from uBAM to fastq using bedtools using the (bamtofastq)[https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html] tool
