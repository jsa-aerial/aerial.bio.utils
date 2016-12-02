;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                 B I O . U T I L S . F I L T E R S                        ;;
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
;; Author: Jon Anthony                                                      ;;
;;                                                                          ;;
;;--------------------------------------------------------------------------;;
;;

(ns aerial.bio.utils.filters

  "Various scoring, filtering, of sequence data. "

  (:require
   [clojure.core.reducers :as r]
   [aerial.utils.coll :refer [vfold] :as coll]
   [aerial.utils.string :as str]
   [aerial.utils.io :refer [letio] :as io]
   [aerial.fs :as fs]
   [aerial.utils.math :as m]
   [aerial.utils.math.probs-stats :as p]
   [aerial.utils.math.infoth :as it]
   [aerial.bio.utils.files :as bufiles]
   [aerial.bio.utils.aligners :as aln]))


(def phred-scores
  "Sangor format Phred score encoding/decoding map. Sangor format
  encodes Phred scores from 0-93 as ASCII using offset 33 (ASCII
  character codes 33-126). NOTE, as of end Feb 2011 all Illumina
  sequencing produces fastqs with this format. "
  (let [pchs "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"]
    (into {} (mapv #(vector %1 %2) pchs (range (count pchs))))))

(defn prob-correct->phred-score
  "Convert a probability p of a base being _correct_ to corresponding
  Phred score. Note that the std mapping (definition of phred score)
  is a mapping from the probability -p of a base being _in_correct."
  [p]
  (- (* 10 (m/log10 (- 1.0 p)))))

(comment
  [(phred-scores \.) (prob-correct->phred-score 0.95)]
  [(phred-scores \.) (prob-correct->phred-score 0.96)])



;;; turn on box warnings
(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)


(defn qcscore-min-entropy
  "Compute minimum entropy for a distribution of kmers of size winsize
  over an msg with 'infinite' non repeated alphabet."
  [^double base%, ^double info%, ^long winsize]
  (let [x (Math/round (* winsize info%))
        infoseq (concat (range (- winsize x)) (range x))]
    [(->> base% (- 1.0) m/log10 double (* 10.0) (-) Math/round)
     ;;(Math/round (- (* 10 (m/log10 (- 1.0 base%)))))
     (->> infoseq (p/probs 1) it/entropy double (* 100.0)
          Math/floor (#(/ (double %) 100)))]))


(defn percent-pass-qscore
  "For qc line qc, return % qc > qc-ctpt value"
  [qc, ^long qc-ctpt]
  (let [len (count qc)
        tfn (comp (map phred-scores) (keep #(when (>= (long %) qc-ctpt) 1)))
        gcnt (reduce (tfn +) 0 qc)]
    (/ (double gcnt) (double len))))


(defn seq-filter
  [[id sq _ qc] & {:keys [qc-ctpt winsize ent-ctpt] :or
      {qc-ctpt 13, winsize 10, ent-ctpt 3.12}}]
  (let [qc-ctpt (long qc-ctpt)
        ent-ctpt (double ent-ctpt)
        cnt (->> qc (map phred-scores)
                 (map #(if (< (long %2) qc-ctpt) \# (char (+ (long %1) 33)))
                      (iterate #(mod (inc (long %)) 90) 0))
                 reverse
                 (coll/sliding-take winsize) (map #(apply str %))
                 (map #(vector % (it/entropy (p/probs 1 %))))
                 (take-while #(< (double (second %)) ent-ctpt)) count)
        gcnt (- (count sq) cnt)
        gsq  (str/take gcnt sq)
        gqc  (str/take gcnt qc)]
    [gcnt id gsq gqc
     (->> gqc (mapv phred-scores)
          (filterv #(>= (long %) qc-ctpt))
          count (#(/ (double %) (double gcnt))))]))


(defn trim-ends
  [[gcnt id ^String gsq gqc qc% :as rec] pretrim min-len sqc% marker]
  (let [gcnt (long gcnt)
        pretrim (long pretrim)
        min-len (long min-len)
        qc% (double qc%)
        sqc% (double sqc%)]
    (if (and (> gcnt min-len) (>= qc% sqc%))
      (let [idx (long (if marker (.indexOf gsq ^String marker min-len) -1))
            gcnt (if (> idx 0) idx gcnt)
            gsq  (str/substring gsq pretrim gcnt)
            gqc  (str/substring gqc pretrim gcnt)
            gcnt (- gcnt pretrim)]
        [gcnt id gsq gqc qc%])
      rec)))


;; Turn off box / reflection warnings
(set! *unchecked-math* false)
(set! *warn-on-reflection* false)


