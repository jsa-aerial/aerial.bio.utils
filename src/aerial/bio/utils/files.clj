;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                   B I O . U T I L S . F I L E S                          ;;
;;                                                                          ;;
;; Permission is hereby granted, free of charge, to any person obtaining    ;;
;; a copy of this software and associated documentation files (the          ;;
;; "Software"), to deal in the Software without restriction, including      ;;
;; without limitation the rights to use, copy, modify, merge, publish,      ;;
;; distribute, sublicense, and/or sell copies of the Software, and to       ;;
;; permit persons to whom the Software is furnished to do so, subject to    ;;
;; the following conditions:                                                ;;
;;                                                                          ;;
;; The above copyright notice and this permission notice shall be           ;;
;; included in all copies or substantial portions of the Software.          ;;
;;                                                                          ;;
;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,          ;;
;; EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF       ;;
;; MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                    ;;
;; NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE   ;;
;; LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION   ;;
;; OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION    ;;
;; WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.          ;;
;;                                                                          ;;
;; Author: Jon Anthony & Shermin Peis                                       ;;
;;                                                                          ;;
;;--------------------------------------------------------------------------;;
;;

(ns aerial.bio.utils.files

  "Various bio sequence file format readers, writers, verifiers, and
   manipulators."

  (:require [clojure.java.io :as cjio]
            [clojure.string :as cstr]
            [clojure.set :as set]
            [clojure.pprint :refer [cl-format]]
            [clojure.data.csv :as csv]

            [aerial.utils.misc :as misc
             :refer [gen-uid getenv raise catch-all with-ckd]]
            [aerial.utils.io :as io
             :refer [letio]]
            [aerial.utils.coll :as coll
             :refer [ensure-vec dropv-until takev in transpose vfold pxmap]]
            [aerial.utils.string :as str]
            [aerial.fs :as fs]
            [aerial.bio.utils.params :as pams]
            [aerial.bio.utils.seqxform
             :refer [reverse-compliment]]
            )
  (:import aerial.utils.io.auioReader))

;;; Tell letio how to work with read-seqs and using auioReader
(aerial.utils.io/add-rdrwtr 'aerial.bio.utils.files/read-seqs true)


;;; ------------------------------------------------------------------------;;;
;;;
;;; Convert Sto and Fasta split sequence format files into conjoined
;;; versions.  Many Sto and Fasta files from various sites come in old
;;; fashioned 80 col mode where sequences are split at 80 column mark.
;;; For aligned files this is even worse as you have groups of
;;; sequences split across lines separated by whole pages of other
;;; sequence (parts).  For example, RFAM alignments.  This group puts
;;; all those back together so that each sequence is on a single line.


(defn write-sto
  "A work in progress...  Write a new sto composed of the various
   given parts to the file spec given as NEWSTO.  AUTH-LINES are the
   authoring header lines - including the STOCKHOLM line.  Generally
   there are two of these - the STOCKHOLM line (with version) and the
   originating author or program that generated the content (for
   example, Infernal).

   COMMENT-LINES is a collection of the #=GF/GC lines, with the
   exception of the GC SS_cons and RF lines.  Comment-lines may be
   empty (for example, []).

   NM-SQ-PAIRS is a collection (typically vector/list) of pairs of the
   entries (name/start-end/strand) and the associated sequence (in
   gapped form).  If this is created via JOIN-STO-FASTA-LINES, the
   vector of [id sq] pairs that is the sequence part of the
   nm-sq-pair, will have the id part filtered out automatically.

   SS-LINES is the set of 'secondary structure' lines.  These are the
   GC SS_cons and RF lines.  SS-LINES may contain the final '//' line
   or not.  If not, it is still written to the file, if so, only the
   one '//' is written.
  "
  [newsto auth-lines comment-lines nm-sq-pairs ss-lines]
  (let [ss-pairs (if (vector? (first ss-lines))
                   (if (vector? (-> ss-lines first second))
                     (map (fn[[gf [id ss]]] [gf ss]) ss-lines)
                     ss-lines)
                   (map #(let [bits (str/split #"\s+" %)]
                           [(cstr/join " " (take 2 bits)) (last bits)])
                        ss-lines))
        nm-sq-pairs (if (vector? (-> nm-sq-pairs first second))
                      (map (fn[[nm [id sq]]] [nm sq]) nm-sq-pairs)
                      nm-sq-pairs)]
    (io/with-out-writer newsto
      (doseq [l auth-lines] (println l))
      (println)
      (doseq [l comment-lines] (println l))
      (when (seq comment-lines) (println))
      (doseq [[nm sq] nm-sq-pairs] (cl-format true "~A~40T~A~%" nm sq))
      (doseq [[gf ss] ss-pairs]
        (if (not= gf "//")
          (cl-format true "~A~40T~A~%" gf ss)
          (println gf)))
      (when (not= (first (last ss-pairs)) "//")
        (println "//")))))


(defn sto-GC-and-seq-lines [stofilespec]
  (coll/separate
   #(and (> (count %) 1)
         (or (not (.startsWith % "#"))
             (or (.startsWith % "#=GC SS_cons")
                 (.startsWith % "#=GC RF"))))
   (filter #(not= (str/replace-re #"\s+" "" %) "")
           (io/read-lines (fs/fullpath stofilespec)))))


(defn join-sto-fasta-lines [infilespec origin]
  (let [[seqcons-lines gc-lines] (sto-GC-and-seq-lines infilespec)
        gc-lines (if (not= origin "")
                   (concat (take 1 gc-lines) [origin] (drop 1 gc-lines))
                   gc-lines)
        recombined-lines (sort-by
                          #(-> % second first)
                          (vec (reduce
                                (fn [m l]
                                  (let [[nm sq]
                                        (cond
                                         ;;splits the line apart and
                                         ;;hopefully creates vector
                                         ;;["#GC SS_cons" structure]
                                         (.startsWith l "#=GC SS_cons")
                                         [(cstr/join " " (butlast (str/split
                                                                  #"\s+" l)))
                                          (last (str/split #"\s+" l))]

                                         (.startsWith l "#")
                                         (str/split #"\s{2,}+" l)

                                         :else
                                         (str/split #"\s+" l))

                                        prev (get m nm [(gen-uid) ""])]
                                    (assoc m  nm [(first prev)
                                                  (str (second prev) sq)])))
                                {} seqcons-lines)))
        {seq-lines false cons-lines true} (group-by
                                           #(or (.startsWith (first %) "//")
                                                (.startsWith (first %) "#"))
                                           recombined-lines)]
    [gc-lines seq-lines cons-lines]))


(defn join-sto-fasta-file
  "Joins (de-blocks) unblocked sequence lines in a sto file or fasta
  file. If in-filespec is a sto file, ORIGIN is a #=GF line indicating
  tool origin of file.  For example, '#=GF AU Infernal 1.0.2'. For
  stos defaults to nothing, for fastas, not used."

  [in-filespec out-filespec
   & {origin :origin :or {origin ""}}]
  (let [fasta? #{"fasta" "fna" "fa" "hitfna"}
        [gc-lines seq-lines cons-lines]
        (join-sto-fasta-lines in-filespec origin)]
    (io/with-out-writer (fs/fullpath out-filespec)
      (if (fasta? (fs/ftype in-filespec))
        (let [sls (->> seq-lines
                       (map first)
                       (partition-by #(re-find #"^>" %))
                       (partition-all 2)
                       (map (fn[[[nm] sbits]] [nm (apply str sbits)])))]
          (doseq [[nm sq] sls]
            (println nm)
            (println (cstr/upper-case sq))))
        ;; Else, sto...
        (do (doseq [gcl gc-lines] (println gcl))
            (doseq [sl seq-lines]
              (let [[nm [_ sq]] sl]
                (cl-format true "~A~40T~A~%" nm sq)))
            (doseq [cl cons-lines]
              (let [[nm [_ sq]] cl]
                (cl-format true "~A~40T~A~%" nm sq))))))))


(defn split-join-fasta-file
  [in-file
   & {:keys [base pat namefn entryfn testfn]
      :or {base "" pat #"^>gi" entryfn identity testfn (fn[x y] true)}}]
  {:pre [(fs/directory? base) (fn? namefn) (fn? testfn)]}
  (doseq [[gi sq] (->> in-file io/read-lines
                       (partition-by #(re-find pat %))
                       (partition-all 2)
                       (map (fn[[[nm] sbits]]
                              [(entryfn nm) (apply str sbits)])))]
    (let [nm (namefn gi)]
      (when (testfn nm sq)
        (io/with-out-writer (fs/join base (str nm ".fna"))
          (println gi)
          (println sq))))))

(defn split-join-ncbi-fasta-file
  "Split a fasta file IN-FILE into the individual sequences and
   unblock the sequence if blocked.  The resulting individual [nm sq]
   pairs are written to files named for the NC name in the gi line of
   in-file and in the DEFAULT-GENOME-FASTA-DIR location.

   The main use of this function is to take a refseq fasta
   db (composed of many multi seq fasta files) and split the db into a
   normed set of named sequence files for quick access to sequence per
   name in various other processing (see gen-name-seq for example).

   Canonical use case example:

   (fs/dodir \"/data2/BioData/Fasta\" ; RefSeqxx fasta files
             #(fs/directory-files % \"fna\")
             #(split-join-ncbi-fasta-file %))
  "
  [in-file]
  (let [base (pams/default-genome-fasta-dir)]
    (doseq [[id sq] (->> in-file io/read-lines
                         (partition-by #(re-find #"^>" %))
                         (partition-all 2)
                         (map (fn[[[nm] sbits]] [nm (apply str sbits)])))]
      (let [nm (->> (str/substring id 1)
                    (str/split #" ") first
                    (str/split #"\.") first)]
        #_(println id nm)
        (when (re-find #"^NC_" nm)
          (io/with-out-writer (fs/join base (str nm ".fna"))
            (println id)
            (println sq)))))))


(defn chunk-genome-fnas
  "Take the set of fnas in directory GENOME-FNA-DIR (presumably
   created by split-join-ncbi-fasta-file or similar) and aggregate
   them into a new smaller set of files, where each new file contains
   the contents of CHUNK-SIZE input files (with the possible exception
   of the last file having a smaller number).  This is useful for
   creating custom data sets for search.
  "
  [genome-fna-dir &
  {:keys [chunk-size] :or {chunk-size 100}}]
  (let [dir (fs/join genome-fna-dir "Chunked")
        all (->> (fs/directory-files genome-fna-dir ".fna") sort
                 (partition-all chunk-size))]
    (when (not (fs/exists? dir)) (fs/mkdir dir))
    (doseq [grp all]
      (let [n1 (-> grp first fs/basename (fs/replace-type ""))
            n2 (-> grp last fs/basename (fs/replace-type ""))
            file (fs/join dir (str n1 "-" n2 ".fna"))]
        (io/with-out-writer file
          (doseq [f grp
                  l (io/read-lines f)]
            (println l)))))))


;;; Convert STO format to ALN format (ClustalW format).  This is
;;; needed by some processors which cannot take a Stockholm alignment
;;; format but need an "equivalent" in ClustalW ALigNment format.
;;;
;;; OK, (9-Feb-2012) some tools seem to need things blocked while
;;; others don't work if they are blocked.  Worse, what counts as
;;; valid Clustal/aln format or not is ill defined with multiple
;;; definitions in the community (e.g., many claim 60 character seqs
;;; per line but others say 50; some claim must be blocked, others say
;;; unblocked is valid).  So, we have two variants.  One which blocks
;;; and a main driver which calls blocked version if requested or just
;;; does simple unblocked itself.
;;;
(defn sto->aln-blocked
  "Convert a stockhom format alignment file into its ClustalW
   equivalent BLOCKED ALN format. Blocking is done in 60 character
   chunks.  STOIN is the filespec for the stockholm format file and
   ALNOUT is the filespec for the resulting conversion (it is
   overwritten if it already exists!)"

  [stoin alnout]
  (let [seq-lines (second (join-sto-fasta-lines stoin ""))
        seq-lines (map (fn [[nm [uid sl]]]
                         [nm [uid (str/partition-stg
                                   60 (str/replace-re #"\." "-" sl))]])
                       seq-lines)]
    (io/with-out-writer alnout
      (println "CLUSTAL W (1.83) multiple sequence alignment\n")
      (loop [x seq-lines]
        (let [[nm [uid sl]] (first x)]
          (when (not-empty sl)
            (do
              (doseq [[nm [uid sl]] x]
                  (cl-format true "~A~40T~A~%" nm (first sl)))
              (println "")
              (recur (map (fn [[nm [uid sl]]]
                            [nm [uid (rest sl)]])
                          x)))))))
    alnout))

(defn sto->aln
  "Convert a stockhom format alignment file into its ClustalW
   equivalent ALN format.  STOIN is the filespec for the stockholm
   format file and ALNOUT is the filespec for the resulting
   conversion (it is overwritten if it already exists!)

   BLOCKED is a boolean indicating whether the output should be
   blocked (60 chars per chunk).  Default is unblocked."

  [stoin alnout & {blocked :blocked :or {blocked false}}]
  (if blocked
    (sto->aln-blocked stoin alnout)
    (let [seq-lines (filter #(not (or (.startsWith % "//") (re-find #"^#" %)))
                            (first (sto-GC-and-seq-lines stoin)))
          seq-lines (map #(str/replace-re #"\." "-" %) seq-lines)]
      (io/with-out-writer (fs/fullpath alnout)
        (println "CLUSTAL W (1.83) multiple sequence alignment\n")
        (doseq [sl seq-lines]
          (println sl)))
      alnout)))




;;; Cool stuff from Shermin, but requires much more refactoring of
;;; various other things of his to make it all work. See his fold-ops.
(defn print-sto
  "takes sequence lines and a structure line and writes it into a sto
  format file. the seq-lines needs to be a collection of [name
  sequence] pairs. structure is a string. Simply prints out to the
  repl."

  [seq-lines structure]
  (println "# STOCKHOLM 1.0\n")
  (doseq [sq seq-lines]
    (let [[nm sq] (if (vector? sq)
                    sq
                    (str/split #"\s+" sq))]
      (cl-format true "~A~40T~A~%" nm (str/replace-re #"\-" "." sq))))
  (cl-format true "~A~40T~A~%" "#=GC SS_cons" structure)
  (println "//"))

#_(defn aln->sto
  "takes an alignment in Clustal W format and produces a sto file by
   using RNAalifold to determine the structure and then making it into
   a sto file adding header and a consensus line"

  [in-aln out-sto & {:keys [fold-alg st]
                     :or {fold-alg "RNAalifold"}}]
  (cond
   (identity st) ;structure provided
   (letio [sq (read-seqs in-aln :type "aln")]
     (io/with-out-writer out-sto
       (print-sto sq st))
     out-sto)

   (= fold-alg "RNAalifold") ;structure from RNAalifold
   (letio [st (fold-aln in-aln)
         sq (read-seqs in-aln :type "aln")]
     (io/with-out-writer out-sto
       (print-sto sq st))
     out-sto) ;return out sto filename

   ;;else use cmfinder
   #_(shell/sh "perl" "/home/kitia/bin/gaisr/src/mod_cmfinder.pl" in-aln out-sto)))



;;; Forward declarations...
(declare
 read-seqs
 entry-parts)

(defn read-fqrec
  "Read a fastq 'record' from a file. IN is an input file
  descriptor (an already opened input-stream reader). Returns a
  quad [id sq aux qc] defining the next fastq record from IN."
  [^java.io.BufferedReader in]
  [(.readLine in)
   (.readLine in)
   (.readLine in)
   (.readLine in)])

(defn read-fqrecs
  "Read n fastq 'records' from a file. IN is an input file
  descriptor (an already opened input-stream reader), and N is the
  number of records (4 line chunks) to read. Returns a vector of
  vector quads [id sq aux qc], each quad representing the id line,
  sequence line, auxilliary line and quality control line (phread
  scores)."
  [in n]
  (loop [recs []
         n n]
    (if (= n 0)
      recs
      (let [rec (read-fqrec in)]
        (if (nil? (rec 0))
          recs
          (recur (conj recs rec)
                 (dec n)))))))

(defn write-fqrec
  "Write a fastq 'record' to a file.  OT is an output file
  descriptor (an already opened output-stream writer).  REC is a
  vector quad [id sq aux qc], representing the id line, the sequence
  line, the auxilliary information line and the quality control line
  for a fastq format file."
  [^java.io.BufferedWriter ot rec]
  (let [[id sq aux qc] rec]
    (.write ot (str id "\n"))
    (.write ot (str sq "\n"))
    (.write ot (str aux "\n"))
    (.write ot (str qc "\n"))))

(defn write-fqrecs
  "Write fastq 'records' to file. OT is an output file descriptor (an
  already opened output-stream writer). RECS is a vector/sequence of
  quads [id sq aux qc], each representing the id line, the sequence
  line, the auxilliary information line and the quality control line
  for a fastq format file."
  [ot recs]
  (doseq [rec recs]
    (write-fqrec ot rec)))


(defn read-farec
  "Read a fasta 'record' from a file. IN is an input file
  descriptor (an already opened input-stream reader). Returns a
  pair [id sq] defining the next fasta record from IN."
  [^java.io.BufferedReader in]
  [(.readLine in)
   (.readLine in)])

(defn read-farecs
  "Read n fasta 'records' from a file. IN is an input file
  descriptor (an already opened input-stream reader), and N is the
  number of records (2 line chunks) to read. Returns a vector of
  vector pairs [id sq], each pair representing the id line and
  sequence line."
  [in n]
  (loop [recs []
         n n]
    (if (= n 0)
      recs
      (let [rec (read-farec in)]
        (if (nil? (rec 0))
          recs
          (recur (conj recs rec)
                 (dec n)))))))
(defn write-farec
  "Write a fasta 'record' to a file.  OT is an output file
  descriptor (an already opened output-stream writer).  REC is a
  vector quad [id sq], representing the id line and the sequence line"
  [^java.io.BufferedWriter ot rec]
  (let [[id sq] rec]
    (.write ot (str id "\n"))
    (.write ot (str sq "\n"))))

(defn write-farecs
  "Write fasta 'records' to file. OT is an output file descriptor (an
  already opened output-stream writer). RECS is a vector/sequence of
  quads [id sq], each representing the id line and the sequence line"
  [ot recs]
  (doseq [rec recs]
    (write-farec ot rec)))


(defn collapse-one
  "Collapse the sequences in fqa, a fastq or fasta file and write the
  collapsed value as a fasta record to fasta, a file spec for fasta
  output."
  [fqa fasta]
  (letio [inf (io/open-file fqa :in)
          otf (io/open-file fasta :out)
          cnt (loop [fqrec (read-fqrec inf)
                     M {}]
                (if (nil? (fqrec 0))
                  (reduce (fn[C [sq cnt]]
                            (let [id (str ">" C "-" cnt)]
                              (write-farec otf [id sq])
                              (inc C)))
                          1 M)
                  (let [sq (fqrec 1)]
                    (recur (read-fqrec inf)
                           (assoc M sq (inc (get M sq 0)))))))]
      [fasta cnt]))

(defn collapse-group
  [pairs]
  (doseq [p pairs]
    (collapse-one (first p) (second p))))


(defn- sample-fxx
  ([p f read-fn write-fn]
   (letio [rdr (io/open-file f :in)]
     (loop [rec (read-fn rdr)
            samps []]
       (if (not (first rec))
         samps
         (if (< (rand) p)
           (recur (read-fn rdr) (conj samps rec))
           (recur (read-fn rdr) samps))))))
  ([p f read-fn write-fn samp-file]
   (letio [rdr (io/open-file f :in)
           wtr (io/open-file samp-file :out)]
     (loop [rec (read-fn rdr)]
       (if (not (first rec))
         samp-file
         (if (< (rand) p)
           (do (write-fn wtr rec)
               (recur (read-fn rdr)))
           (recur (read-fn rdr))))))))

(defn sample-fna
  "Sample the sequences in f, a fasta file, with probability
  p. Returns a seq of pairs suitable for writing a fasta file: [id,
  sq]. The id is the corresponding id of the sq in f. In the 3 arg
  case, sampfna is a filespec for an output fasta file where the
  sampling is written."
  ([p f]
   (sample-fxx p f read-farec write-farec))
  ([p f sampfa]
   (sample-fxx p f read-farec write-farec sampfa)))

(defn sample-fq
  "Sample the sequences in f, a fastq file, with probability
  p. Returns a seq of quadtuples suitable for writing a fastq
  file: [id, sq, qcdesc qc]. The id is the corresponding id of the sq
  in f. qcdesc and qc are the corresponding quality description line
  and the quality score line. In the 3 arg case, sampfq is a filespec
  for an output fastq file where the sampling is written."
  ([p f]
   (sample-fxx p f read-fqrec write-fqrec))
  ([p f sampfq]
   (sample-fxx p f read-fqrec write-fqrec sampfq)))

(defn rxy
  "Obtain the `lane-cell-x-y' coordinates of an Illumina fastq read.
  These are obtained from the header information defined for both
  bcl2fastq and bcl-convert output. `hdrl' is the header line of an
  Illumina output fastq read record."
  [hdrl]
  (let [rv (str/split hdrl #":")
        [lane cell x] (->> rv (drop 3) (take 3))
        y (-> 6 rv (str/split #" ") first)]
    (str lane "-" cell "-" x "-" y)))

(defn sampR2fq
  "Sample the reads in the read 2 (R2) fastq file `R2-infq' that is the
  pair of the original read 1 (R1) fastq file that produced the
  sampling in the `R1sampfq' fastq file. The reads in R2-infq
  corresponding to those in R1sampfq are written to `R2-otfq'."
  [R1sampfq R2-infq R2-otfq]
  (letio [R1in (io/open-file R1sampfq :in)
          R2in (io/open-file R2-infq :in)
          R2ot (io/open-file R2-otfq :out)]
    (loop [R1rec (read-fqrec R1in)
           cnt 0]
      (if (nil? (R1rec 0))
        cnt
        (let [R1xy (rxy (R1rec 0))]
          (loop [R2rec (read-fqrec R2in)]
            (if (= R1xy (rxy (R2rec 0)))
              (write-fqrec R2ot R2rec)
              (recur (read-fqrec R2in))))
          (recur (read-fqrec R1in)
                 (inc cnt)))))))

(defn sample-paired-fqs
  "Sample the **paired** reads in R1fq and R2fq with probability
  `p'. The pairing of the reads is maintained, and the output
  files (for R1fq and R2fq) are written to output directory `otdir'."
  [p [R1fq R2fq] otdir]
  (let [R1otfq (->> R1fq fs/basename (fs/join (fs/fullpath otdir)))
        R2otfq (->> R2fq fs/basename (fs/join (fs/fullpath otdir)))]
    (sample-fq p R1fq R1otfq)
    (sampR2fq R1otfq R2fq R2otfq)))


(defn fastq->fna
  "Convert a fastq format file to a fasta file. Fastq files have 4
  lines per 'record' (id, sq, qcid, qc/phred-scores), while the
  corresponding fasta has only the id and sq lines per 'record'. Both
  FQ and FAOT are file specs."
  [fq faot]
  (letio [rdr (io/open-file fq :in)
          wtr (io/open-file faot :out)]
    (loop [rec (read-fqrec rdr)]
      (if (not (rec 0))
        faot
        (do (write-farec wtr [(str ">" (rec 0)) (rec 1)])
            (recur (read-fqrec rdr)))))))

(defn fastqs->fnas
  "Convert the fastqs in FQS (a seq) to corresponding fasta files. The
  fasta files are named as the fastqs but with file type .fna. If
  outdir is given, place the fastas there, otherwise place in same
  directory as corresponding fastqs. The items in fqs are file specs!"
  [fqs & {:keys [outdir]}]
  (let [outdir (if outdir outdir (fs/dirname (or (first fqs) "")))]
    (doseq [fq fqs]
      (let [fa (fs/join outdir (-> fq fs/basename (fs/replace-type ".fna")))]
        (fastq->fna fq fa)))))


(defn fna->fastq
  "Convert a fasta format file into a fastq format file. infna and
  outfq are fie specifications (strins) of the input fasta and the
  output fastq.

  The issue here is simply what should the quality control (Phred
  scores) be? We assume the user understands the issue and provide a
  'slight' helper for this: bq% is the probability the bases (all!!)
  are correct. We don't accommodate the case of providing bqs for each
  base, as that would indicate you already have a fastq version of the
  sq in question (or certainly _should_ have such). The bc% is
  converted to the corresponding Phred score and then Sanger
  encoded (0 -> 33, and all such scores written as their corresponding
  ASCII character.)

  NOTE: output is unblocked even if input fasta is in blocked format!

  "
  [infna outfq bc%]
  (io/with-out-writer (fs/fullpath outfq)
    (let [[_ seq-lines] (join-sto-fasta-lines infna "")
          nmsq-pairs (->> seq-lines
                          (map first)
                          (partition-by #(re-find #"^>" %))
                          (partition-all 2)
                          (map (fn[[[nm] sqbits]] [nm (apply str sqbits)])))
          phred (- (* 10 (aerial.utils.math/log10 (- 1.0 bc%))))
          sanger (-> phred (+ 33) long char)]
      (doseq [[nm sq] nmsq-pairs]
        (println (str "@" (str/substring nm 1)))
        (println (cstr/upper-case sq))
        (println (str "+ Generated by fasta->fastq with QC scores for " bc%
                      " probability base good"))
        (dotimes [_ (count sq)] (print sanger))
        (println)))))

(defn fnas->fastqs
  "Convert fastas (collection of fasta filespecs as strings) to
  corresponding fastqs in directory fastqdir. Conversion is by
  fna->fastq."
  [fastas fastqdir bc%]
  (doseq [fa fastas]
    (let [fq (fs/join fastqdir (-> fa fs/basename (fs/replace-type ".fq")))]
      (fna->fastq fa fq bc%))))


(defn gbank->fna
  "Obtains the genomic sequence (the 'ORIGIN' record) of a genbank
  format file gbfile and writes it as the sequence (unblocked) of a
  fasta file located in fafile. The id line of the fasta consists of
  the 'LOCUS' fields of the gbfile."
  [gbfile fafile]
  (letio [lines (io/read-lines gbfile)]
    (let [ncid (->> lines first (str/split #"\s+")
                    rest (clojure.string/join " ")
                    (str ">"))
          lines (rest (coll/drop-until #(re-find #"^ORIGIN" %) lines))]
      (loop [lines lines
             sq []]
        (if (seq lines)
          (recur (rest lines)
                 (conj sq (->> lines first (str/split #"\s+")
                               (coll/dropv 2)
                               (apply str) clojure.string/upper-case)))
          (io/with-out-writer fafile
            (println ncid)
            (println (clojure.string/join sq))))))))

(defn gbanks->fnas
  "Calls gbank->fna on all genbank format files in gbfiles, a
  collection of fully qualified file specifications of genbank
  files. The corresponding fasta file specification is composed of the
  basename of a genbank file, with '.fna' file type and placed in
  outdir. Outdir is a string which is an output directory
  specification, which must exist and be a directory."
  [gbfiles outdir]
  (assert (and (string? outdir) (fs/exists? outdir) (fs/directory? outdir))
          (str "Non existent directory " outdir ))
  (doseq [gbf gbfiles]
    (let [fa (fs/join outdir (-> gbf fs/basename (fs/replace-type ".fna")))]
      (gbank->fna gbf fa))))


(defn gbank-loci-sane-loci
  "GBFF loci are often, if not incoherent, certainly vague and not very
  useful. This function takes such and turns them into nice clean loci
  with simple start, end and strand. BUG: this will lose some
  information in certain cases, but should not with CDS and genes. "
  [gbloci]
  (let [std (if (re-find #"complement" gbloci) "-1" "1")]
    (->> gbloci (str/replace-re #"(<|>)" "")
         (re-seq #"[0-9]+\.\.[0-9]+")
         (mapv #(conj (str/split #"\.+" %) std)))))

(defn get-full-comma-sep-stg [stg rdr]
  (loop [l (.readLine rdr)
         stg stg]
    (let [newstg (->> l clojure.string/trim (str stg))]
      (if (not= (str/get l (-> l count dec)) \,)
        newstg
        (recur (.readLine rdr) newstg)))))

(defn genbank-recs
  "Takes a genbank file, gbfile - a string file specification, and a set
  of features (feats) and associated attributes (attrs), so called
  qualifiers and values, and returns a vector where the first element
  describes the LOCUS of the gbfile and each subsequent element is a
  triple [feature loci attr-map], where feature (a string) is one of
  feats, loci is a vector [start end strand], (all strings, strand =
  1|-1), and attr-map is a map with keys in attrs and their associated
  values from the genbank record."
  [gbfile & {:keys [feats attrs]
             :or {feats ["gene", "CDS", "misc_feature",
                         "rRNA", "tRNA", "tmRNA", "ncRNA"]
                  attrs ["gene", "locus_tag", "old_locus_tag", "protein_id",
                         "condon_start", "gene_synonym",
                         "db_xref"]}}]
  (letio [feats (set feats)
          attrs (set attrs)
          trim clojure.string/trim
          blank? clojure.string/blank?
          feature? (fn[l] (and (not (blank? (str/substring l 0 6)))
                              (feats (->> l trim (str/split #"\s+") first))))
          rdr (clojure.java.io/reader gbfile)
          locus (->> rdr .readLine (str/split #"\s+") vec)]
    (while (->> rdr .readLine (re-find #"^FEATURE") not))
    (loop [l (.readLine rdr)
           recs [locus]
           cnt 0]
      (cond
        (not= (str/get l 0) \space)
        recs

        (feature? l)
        (let [[feat loci] (->> l trim (str/split #"\s+"))
              loci (if (= (str/get loci (-> loci count dec)) \,)
                     (get-full-comma-sep-stg loci rdr)
                     loci)
              loci (gbank-loci-sane-loci loci)
              [rec l]
              (loop [l (.readLine rdr)
                     fields {}]
                (if (not (blank? (str/substring l 0 6)))
                  [[feat loci fields] l]
                  (let [[field val] (->> l trim (str/split #"="))
                        field (when (not= field "") (str/substring field 1))]
                    (if (attrs field)
                      (recur (.readLine rdr)
                             (assoc fields field val))
                      (recur (.readLine rdr)
                             fields)))))]
          (recur l (conj recs rec) (inc cnt)))

        :else
        (recur (.readLine rdr)
               recs
               cnt)))))

(defn genbk2gtf
  "Extract a GTF (gene transfer format) file from the record
  information (features and their attributes) given in the genbank
  format file gbfile. The GTF file also includes the p_id attribute on
  CDS records, as required by the cuff* suite of software (this is no
  longer all that meaningful as the cuff* suite of tools is basically
  deprecated as 'not good enough').

  options is a map of options to use. Currently supports
  keys :id-order and :feats.

  If :id-order is given, the value should be a vector of
  \"locus_tag\", \"old_locus_tag\", \"gene\", and \"protein_id\", the
  order given will determine which of these is used for gene_id and
  transcript_id. The default value is [\"locus_tag\",
  \"old_locus_tag\", \"gene\", and \"protein_id\"]

  If :feats is given, the value should be a vector of feature type
  names taken from this list:

  \"gene\", \"CDS\", \"misc_feature\", \"rRNA\", \"tRNA\", \"tmRNA\", \"ncRNA\"

  This will determine what features to include and, importantly, the
  order will determine the feature type to encode them as. Features in
  genbank files are often encoded under multiple types, for example,
  rRNAs are also listed as genes. If you gave :feats [\"CDS\" \"rRNA\"
  \"gene\"] the rRNA records will have a feature type of rRNA instead
  of gene as rRNA occurs before gene in the list.
  "
  ([gbfile gtfout] (genbk2gtf gbfile gtfout {}))
  ([gbfile gtfout options]
   (let [format clojure.pprint/cl-format
         {:keys [id-order feats]
          :or {id-order ["locus_tag" "old_locus_tag" "gene" "protein_id"]
               feats ["CDS" "rRNA" "gene"]}} options
         recs (genbank-recs gbfile :feats feats)
         acc (->> recs first second)
         src "aerial"
         loci-map (reduce
                   (fn[M [k V m]]
                     (reduce (fn [M v]
                               (assoc M v (assoc (get M v {}) k m)))
                             M V))
                   {} (rest recs))
         gtfrecs (->> loci-map (sort-by (fn[[[s e st] x]] (Integer. s)))
                      (map (fn[[k v]]
                             [k (reduce (fn[x f]
                                          (if (v f) (reduced [f (v f)]) x))
                                        :na feats)]))
                      (map (fn[[[s e st] v]] [(first v) s e st (second v)])))
         p_id (volatile! 0)]
     (io/with-out-writer gtfout
       (doseq [[feat s e st m] gtfrecs]
         (let [id (->> id-order (map #(m %)) (coll/dropv-until identity) first)
               pid (if (= feat "CDS") (str "\"P" (vswap! p_id #(inc %)) \") nil)
               gene (m "gene")
               ptn (m "protein_id")
               attrs (format nil "gene_id ~A; transcript_id ~A" id id)
               ;;attrs (if gene (format nil "~A; gene ~A" attrs gene) attrs)
               attrs (if ptn (format nil "~A; protein_id ~A" attrs ptn) attrs)
               attrs (if pid (format nil "~A; p_id ~A" attrs pid) attrs)
               strand (if (= st "1") "+" "-")]
           (println (format nil "~A\t~A\t~A\t~A\t~A\t.\t~A\t0\t~A"
                            acc src feat s e strand attrs))))))))


(defn sto->fna
  "Convert a sto file into a fasta file.  Split seq lines into names
   and seq data and interleave these.  Seq data has all gap characters
   removed."
  [stoin fnaout]
  (letio []
    (io/with-out-writer fnaout
      (doseq [[nm seqdata] (read-seqs stoin :info :both)]
        (println (str ">" nm))
        (println (str/replace-re #"[-.]+" "" seqdata)))))
  fnaout)


(defn check-sto
  "Checks a sto file to ensure that there are valid characters being
   used in the sequences consensus structure line. Will print out
   errors in the sto file by sequence number.  Input requires a sto
   file"

  [sto & {printem :printem :or {printem true}}]
  (let [valid-symbols #{"A" "C" "G" "U"
                        "-" "." ":" "_" ","
                        "a" "b" "B" "n" "N"
                        "(" ")" "<" ">"}
        [_ seq-lines cons-lines] (join-sto-fasta-lines sto "")
        [_ [_ cl]] (first cons-lines)
        ;;adds namez to the sequences so that they
        ;;can be identified
        sl (map (fn [[nm [_ sq]]]
                  [nm (.toUpperCase sq)])
                seq-lines)

        ;;checks for valid symbols
        check-char (fn [[n s]]
                     [n (every? #(contains? valid-symbols %)
                                (rest (str/split #"" s)))])
        ;;check for common case of repeat named seqs
        check-double-len (fn [[n s]]
                           [n (>= (count s) (* 2 (count cl)))])
        ;;checks to see if sequences have same
        ;;length as consensus structure
        check-len (fn [[n s]]
                    [n (= (count s) (count cl))])

        ;; Check for valid entry formats.  Defined as parseable by
        ;; entry-parts
        check-entry (fn [[n _]]
                      [n (not (vector? (-> n entry-parts catch-all)))])

        chks [["sequence contains invalid character in: "
               (map first (remove #(second %) (map #(check-char %) sl)))]
              ["sequence is repeated - two or more times with same name: "
               (map first (filter #(second %) (map #(check-double-len %) sl)))]
              ["sequence contains invalid length compared to cons-line in: "
               (map first (remove #(second %) (map #(check-len %) sl)))]
              ["sequence line has invalid entry - not parseable: "
               (map first (remove #(not (second %)) (map #(check-entry %) sl)))]
              ["consensus structure contains invalid character"
               (if-let [x (second (check-char ["#SS_cons" cl]))]
                 ()
                 ["#SS_cons"])]]]

    (when printem
      (doseq [[msg x] chks]
        (when (seq x)
          (prn msg x))))

    (if (every? #(-> % second seq not) chks)
      (do (when printem (prn "sto is good"))
          :good)
      (filter #(seq (second %)) chks))))


;;; ----------------------------------------------------------------------

(defn make-entry
  "'Inverse' of entry-parts.  EVEC is a vector of shape [nm [s e]
   strd], where nm is the entry name, S and E are the start and end
   coordinates in the genome, and strd is the strand marker, 1 or
   -1. Returns the full entry as: nm/s-e/strd
  "
  ([evec]
     (let [[nm [s e] st] evec]
       (str nm "/" s "-" e "/" st)))
  ([nm s e st]
     (make-entry [nm [s e] st])))


(defn entry-parts
  "ENTRY is a string \"name/range/strand\", where name is a genome
   name, range is of the form start-end and strand is 1 or -1.  At
   least name must be supplied.  DELTA is an integer, which will be
   subtracted from start and added to end.

   Returns a triple [name [start end] strand]
  "
  [entry & {:keys [ldelta rdelta] :or {ldelta 0 rdelta 0}}]
  (let [[name range] (str/split #"( |/|:)+" 2 entry)
        name (re-find #"[A-Za-z0-9]+_[A-Za-z0-9_-]+" name)
        [range strand] (if range (str/split #"/" range) [nil nil])
        [s e st] (if (not range)
                   [1 Long/MAX_VALUE "1"]
                   (let [range (if (= \c (.charAt range 0))
                                 (subs range 1)
                                 range)
                         [s e] (map #(Integer. %) (str/split #"-" range))
                         strand (if strand strand (if (> s e) "-1" "1"))
                         [s e] (if (< s e) [s e] [e s])
                         [s e] [(- s ldelta) (+ e rdelta)]
                         [s e] [(if (<= s 0) 1 s) e]]
                     [s e strand]))]
    [name [s e] st]))


(defn gen-entry-file
  "Take entries (a collection of nm/s-e/std - see make-entry /
  entry-parts) and write them out to file. Returns file"
  [entries file]
  (io/with-out-writer (fs/fullpath file)
    (doseq [e entries]
      (println e)))
  file)

(defn gen-entry-nv-file
  "Take entries (a collection of entry [,v]* where v(s) are optional
  values (JSD, EVals, etc)) and write them as csv file."
  [entries file]
  (io/with-out-writer (fs/fullpath file)
    (doseq [e entries]
      (if (coll? e)
        (println (->> e (map str) (cstr/join ", ")))
        (println e))))
  file)


(defn has-loc? [entries]
  "Return whether entries elements have loci (full entries)"
  (let [entries (cond
                 (seq? entries) entries
                 (fs/exists? entries) (str/split #"\n" (slurp entries))
                 :else (raise :type :unknown-entries :entries entries))]
    (let [e (first (ensure-vec entries))]
      (re-find #"[0-9]+\-[0-9]+" e))))




(defn gen-name-seq
  "Generate a pair [entry genome-seq], from ENTRY as possibly modified
   by [L|R]DELTA and RNA.  ENTRY is a string \"name/range/strand\",
   where

   name is the genome NC name (and we only currently support NCs),

   range is of the form start-end, where start and end are integers (1
   based...) for the start and end coordinates of name's sequence to
   return.  NOTE: start < end, as reverse compliment information comes
   from strand.

   strand is either -1 for reverse compliment or 1 for standard 5'->3'

   LDELTA is an integer, which will be _subtracted_ from start.  So,
   ldelta < 0 _removes_ |ldelta| bases from 5', ldelta > 0 'tacks on'
   ldelta extra bases to the 5' end.  Defaults to 0 (no change).

   RDELTA is an integer, which will be _added_ to the end.  So, rdelta
   < 0 _removes_ |rdelta| bases from 3', rdelta > 0 'tacks on' rdelta
   extra bases to the 3' end.  Defaults to 0 (no change).

   If RNA is true, change Ts to Us, otherwise return unmodified
   sequence.

   BASEDIR is the location of the NC fasta files.  Generally, this
   should always be the default location.
  "
  [entry & {:keys [basedir ldelta rdelta rna]
            :or {basedir (pams/default-genome-fasta-dir entry)
                 ldelta 0 rdelta 0 rna true}}]
  (letio [basedir (if (fn? basedir) (basedir entry) basedir)
          [name [s e] strand] (entry-parts entry :ldelta ldelta :rdelta rdelta)
          _ (when (<= s 0)
              (raise :type :bad-entry
                     :input [entry ldelta rdelta]
                     :result [name s e strand]))
          entry (str name "/" s "-" e "/" strand)
          fname (fs/join basedir (str name ".fna"))
          sq (->> (io/read-lines fname) second
                  (str/drop (dec s)) (str/take (- (inc e) s)))
          sq (if (= strand "-1") (reverse-compliment sq) sq)
          sq (if rna (str/replace-re #"T" "U" sq) sq)]
    [entry sq]))


;;; (fs/dodir
;;;  "/home/kaila/Bio/Work/S15work/062612"
;;;  #(fs/directory-files % ".ent")
;;;  gen-name-seq-pairs)
;;;
(defn gen-name-seq-pairs
  "Generate a sequence of pairs [name genome-seq] from the given
   collection denoted by entries, which is either a collection or a
   string denoting an entry filespec (either hand made or via
   gen-entry-file or similar).

   See gen-name-seq for details.  Basically this is (map gen-name-seq
   entries) with some options.  In particular, entries may lack a
   strand component (again see gen-name-seq for format details), and
   STRAND here would be used to indicate all entries are to have this
   strand.  If entries have a strand component, strand should be 0,
   otherwise you will force all entries to the strand given.  Strand
   is either -1 (reverse compliment) or 1 for standard 5'->3', or 0
   meaning 'use strand component in entries'

   BASEDIR is the location of the NC fasta files.  Generally, this
   should always be the default location.
  "
  [entries & {:keys [basedir strand ldelta rdelta rna]
              :or {basedir pams/default-genome-fasta-dir
                   strand 0 ldelta 0 rdelta 0 rna true}}]
  (let [entries (if (string? entries) (io/read-lines entries) entries)
        entries (if (= 0 strand) entries (map #(str % "/" strand) entries))]
    (map #(gen-name-seq
           % :basedir basedir :ldelta ldelta :rdelta rdelta :rna rna)
         entries)))




(defn nms-sqs->fasta-file [nms-sqs filespec]
  (let [filespec (fs/fullpath filespec)]
    (io/with-out-writer filespec
      (doseq [[nm sq] nms-sqs]
        (println (str ">" nm))
        (println sq)))
    filespec))


(defn entry-file->fasta-file
  [efile & {:keys [names-only]}]
  (letio [efile (fs/fullpath efile)
          fasta-filespec (fs/fullpath (fs/replace-type efile ".fna"))
          entseqs (if names-only
                    (read-seqs efile :info :both)
                    (gen-name-seq-pairs (io/read-lines efile)))]
    (io/with-out-writer fasta-filespec
      (doseq [[entry sq] entseqs]
        (println (str ">" entry))
        (println sq)))
    fasta-filespec))


(defn get-selection-fna
  "Get a fasta file for SELECTIONS.  If selections is a file, assumes
   it is in fact the fasta file.  If selections is a collection,
   assumes the collection is a set of pairs [nm sq], and converts to
   corresponding fasta file.  Returns full path of result file."
  [selections]
  (if (coll? selections)
    (nms-sqs->fasta-file selections)
    (fs/fullpath selections)))


(defn gen-nc-genome-fnas
  "One shot genome sequence fnas generator.  Typically used once per
   data update.  Needs to be generalized to be able to use Genbank fna
   archives.  Currently assumes one fna AND ONLY NC_* genomes.
  "
  [full-nc-genomes-filespec]
  (let [dir (fs/dirname full-nc-genomes-filespec)
        fna-pairs (partition-all 2 (io/read-lines full-nc-genomes-filespec))]
    (doseq [[gi genome] fna-pairs]
      (let [fname (fs/join dir (str (re-find #"NC_[0-9]+" gi) ".fna"))]
        (io/with-out-writer fname
          (println gi)
          (println genome))))))




;;; ----------------------------------------------------------------------- ;;;

;;;  This crap needs to be refactored and most of it probably
;;;  eliminated.  Most of it is old stuff that is largely superceded
;;;  but still required by things in post-db-csv, but we at least fix
;;;  up the names of the csv processors to more accurately reflect
;;;  what they are!


(defn get-legacy-csv-info [rows]
  (let [newnes [:rfam :overlap :new]]
    (loop [entries []
           rows (drop 1 rows)]
      (let [entry (first rows)]
        (if (< (count entry) 24)
          (sort #(str/string-less? (%1 0) (%2 0)) entries)
          (recur (conj entries
                       [(entry 4) (entry 9) (entry 10) (entry 24)
                        (entry 11)
                        (newnes (Integer/parseInt (entry 15)))
                        (entry 13)])
                 (drop 1 rows)))))))

(defn get-gaisr-csv-info [rows]
  (loop [entries []
         rows (drop 1 rows)]
    (let [entry (first rows)]
      (if (< (count entry) 12) ; minimum used fields
        (sort #(str/string-less? (%1 0) (%2 0)) entries)
        (recur (conj entries
                     [(entry 0) (entry 3) (entry 4) (entry 9)
                      (entry 8)
                      :new "N/A"])
               (drop 1 rows))))))


(defn canonical-csv-entry-info [entries & {:keys [ev] :or {ev 0.0}}]
  (map #(let [[nm [s e] sd] (entry-parts %)
              [s e] (if (= sd "1") [s e] [e s])]
          [nm s e ev 0.0 :new sd])
       entries))

(defn get-sto-as-csv-info [stofile & {:keys [ev] :or {ev 0.0}}]
  (let [entries (read-seqs stofile :info :name)]
    (canonical-csv-entry-info entries :ev ev)))

(defn get-ent-as-csv-info [ent-file]
  (let [entries (io/read-lines ent-file)]
    (if (> (->> entries first csv/read-csv first count) 1)
      (let [einfo (canonical-csv-entry-info
                   (map #(->> % (str/split #",") first) entries))
            entropy-scores (->> ent-file slurp csv/read-csv
                                butlast (map second))]
        (map (fn[[nm s e _ _ x y] score] [nm s e score 0.0 x y])
             einfo entropy-scores))
      (canonical-csv-entry-info entries))))


(defn get-csv-entry-info [csv-hit-file]
  (let [file (fs/fullpath csv-hit-file)
        ftype (fs/ftype file)]
    (cond
     (= ftype "sto") (get-sto-as-csv-info file)
     (= ftype "ent") (get-ent-as-csv-info file)
     :else
     (let [rows (csv/read-csv (slurp file))
           head (first rows)]
       (if (= (first head) "gaisr name")
         (get-gaisr-csv-info rows)
         (get-legacy-csv-info rows))))))



;;; This one was from dists and punted to the old get-enties in
;;; post-db-csv for csv entries
;;;
(defn get-entries
  [filespec & [seqs]]
  (let [fspec (fs/fullpath filespec)
        ftype (fs/ftype fspec)]
    (if (not= ftype "csv")
      (read-seqs filespec :info (if seqs :data :name))
      (->> (get-csv-entry-info fspec)
           (keep (fn[[nm s e sd]]
                   (when (fs/exists? (fs/join (pams/default-genome-fasta-dir nm)
                                              (str nm ".fna")))
                     (str nm "/"
                          (if (> (Integer. s) (Integer. e))
                            (str e "-" s "/-1")
                            (str s "-" e "/1"))))))
           (#(if seqs (map second (gen-name-seq-pairs %)) %))))))




;;;------------------------------------------------------------------------;;;
;;; Various sequence file readers for various formats.  In particular,
;;; fna, sto, and aln


(defn gaisr-seq-set?
  "Returns true if either

   X is a filespec string with type extension fna, aln, sto, or gma
   X is a java.io.File (presumed to be of one of the above formats)
   X is a collection (presumed to have seqences as elements)
  "
  [x]
  (let [ftypes #{"fna" "fa" "hitfna" "aln" "sto" "gma"}]
    (or (and (string? x)
             (ftypes (fs/ftype x)))
        (isa? (type x) java.io.File)
        (coll? x))))


(defn seqline-info-mapper
  "Helper function for READ-SEQS.  Returns the function to map over
   seq lines to obtain the requested info.  TYPE is supported seq file
   type (aln, sto, fna, fa, gma).  INFO is either :name for the
   sequence identifier, :data for the sequence data, or :both for name
   and data.

   Impl Note: while this almost begs for multimethods, that would
   actually increase the complexity as it would mean 14 methods to
   cover the cases...
  "
  [type info]
  (if (= info :both)
    #(do [((seqline-info-mapper type :name) %)
          ((seqline-info-mapper type :data) %)])
    (case type
          "aln"
          (if (= info :data)
            #(str/replace-re #"^N[CZ_0-9]+\s+" "" %)
            #(re-find  #"^N[CZ_0-9]+" %))

          "sto"
          (if (= info :data)
            #(str/replace-re #"(N[CZ_0-9]+|[A-Za-z0-9._/-]+)[,\s]+" "" %)
            #(second (re-find  #"(N[CZ_0-9]+|[A-Za-z0-9._/-]+)[,\s]+" %)))

          "ent"
          ;; This is a bit annoying as we need to account for ent
          ;; files with csv annotation beyond entries
          (if (= info :data)
            #(-> % csv/read-csv ffirst gen-name-seq second)
            #(->> (str/split #"([\s,]|/)" %) (take 3) (cstr/join "/")))

          "gma" (raise :type :NYI :info "GMA format not yet implemented")

          ("fna" "fa" "hitfna" "fasta")
          (if (= info :data)
            second
            #(str/replace-re #">+" "" (first %))))))

(defn read-seqs
  "Read the sequences in FILESPEC and return set as a lazy seq.
   Filespec can denote either a fna, fa, hitfna, aln, sto, or gma file
   format file.
  "
  [input & {:keys [info ftype] :or {info :data}}]
  #_(prn input)
  (let [[file rdr] (cond
                     (= (type input) auioReader)
                     [(:fspec input) (:rdr input)]
                     ;; Should be fspec, which works as rdr for read-lines
                     (string? input)
                     [input input]
                     ;; Sould be reader and must have type given
                     :else
                     [nil input])
          ftype (or ftype (fs/ftype file))
          lines (io/read-lines rdr)
          sqs (filter #(let [l (str/replace-re #"^\s+" "" %)]
                         (and (not= l "") (not (.startsWith l "#"))))
                      lines)

          fun (seqline-info-mapper ftype info)
          sqs (if (re-find #"^CLUSTAL" (first sqs)) (rest sqs) sqs)
          sqs (dropv-until #(re-find #"^(>|)[A-Za-z0-9@]" %) sqs)
          sqs (if (#{"fna" "fa" "hitfna" "fasta"} ftype) (partition 2 sqs) sqs)
        sqs (if (= ftype "sto")
              (take-while #(re-find #"^[0-9A-Z]" %) sqs)
              sqs)]
    (map fun sqs)))


(defn read-aln-seqs
  "Read the _aligned_ sequences in FILESPEC and return them in a
  Clojure seq.  Filespec can denote either an aln, sto, or gma file
  format file.  If COLS is true, return the _columns_ of the
  alignment (including gap characters).
  "
  [filespec & {cols :cols :or {cols false}}]
  {:pre [(in (fs/ftype filespec) ["aln" "sto" "gma"])]}
  (let [sqs (read-seqs filespec)]
    (if cols (transpose sqs) sqs)))

(defn read-dirs-aln-seqs
  "Apply read-aln-seqs across FTYPES in DIR.  FTYPES is a vector of
   one or more file formats (as type extensions) taken from #{\"fna\",
   \"aln\", \"sto\", \"gma\"}, producing a seq of sequence sets, each
   set being a seq of the sequences in a matching file.  NOTE: each
   such set can be viewed as a 'matrix' of the elements (bases/gaps)
   of the sequences.

   If COLS is true, return the transpose of each obtained sequence
   'matrix', i.e., return the seq of column seqs.  NOTE: as the
   formats are alignment formats, all sequences in a file are of the
   same length (including gaps).

   If DIRDIR is true, dir is taken as a directory of directories, each
   of which will have read-aln-seqs applied to it per the above
   description and the return will be a seq of all such applications.

   So, if dirdir is false, the result will have the form:

   (seqs-from-file1 ... seqs-from-filen), where filei is in dir

   If dirdir is true, the result will be nested one level more:

   ((seqs-from-dir1-file1 ... seqs-from-dir1-filek)
    ...
    (seqs-from-dirn-file1 ... seqs-from-dirn-filel))
  "
  [dir & {dirdir :dirdir cols :cols ftypes :ftypes
          :or {dirdir false cols false ftypes ["sto"]}}]
  (let [one-dir (fn[dir]
                  (filter #(> (count %) 1)
                          (fs/dodir
                           dir (fn[d]
                                 (apply set/union
                                        (map #(fs/directory-files d %)
                                             ftypes)))
                           read-aln-seqs)))
        sqs (if dirdir
              (map one-dir (sort (fs/directory-files dir "")))
              (one-dir dir))]
    (if cols
      (if dirdir
        (map #(map transpose %) sqs)
        (map transpose sqs))
      sqs)))


(defn map-aln-seqs
  ""
  ([f cols filespec]
     (list (f (read-aln-seqs filespec :cols cols))))
  ([f par cols filespec & filespecs]
     (pxmap f par (map #(read-aln-seqs % :cols cols)
                       (cons filespec filespecs)))))

(defn reduce-aln-seqs
  ""
  ([f fr cols filespecs]
     (let [par (if (and (coll? filespecs) (> (count filespecs) 4)) 4 1)]
       (reduce fr (apply map-aln-seqs f par cols filespecs))))
  ([f fr v cols filespecs]
     (let [par (if (and (coll? filespecs) (> (count filespecs) 4)) 4 1)]
       (reduce fr v (apply map-aln-seqs f par cols filespecs)))))


(defn map-seqs
  ([f filespec]
   (list (f (read-seqs filespec))))
  ([f par filespec & filespecs]
   ;;(prn :par par :filespec filespec :filespecs filespecs)
   (pxmap f par (map #(read-seqs %)
                     (cons filespec filespecs)))))

(defn reduce-seqs
  ""
  ([f fr filespecs]
     (let [par (if (and (coll? filespecs) (> (count filespecs) 4)) 4 1)]
       (reduce fr (apply map-seqs f par filespecs))))
  ([f fr v filespecs]
     (let [par (if (and (coll? filespecs) (> (count filespecs) 4)) 4 1)]
       (reduce fr v (apply map-seqs f par filespecs)))))


