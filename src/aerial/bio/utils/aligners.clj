;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                  B I O . U T I L S . A L I G N E R S                     ;;
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

(ns aerial.bio.utils.aligners

  "Various sequence aligners by means of low level hackery..."

  [:require
   [clojure.core.reducers :as r]
   [clojure.string :as cljstr]
   [aerial.utils.misc :as aum]
   [aerial.utils.string :as str]
   [aerial.utils.coll :as coll]])


(def nt4submat (into {} (for [c1 [\a \c \g \t]
                              c2 [\a \c \g \t]]
                          [[c1 c2] (if (= c1 c2) 5 -2)])))

(defmacro deep-aget
  ([hint array idx]
    `(aget ~(vary-meta array assoc :tag hint) ~idx))
  ([hint array idx & idxs]
    `(let [a# (aget ~(vary-meta array assoc :tag 'objects) ~idx)]
       (deep-aget ~hint a# ~@idxs))))

(defmacro deep-aset [hint array & idxsv]
  (let [hints '{doubles double, longs long, ints int}
                ;;vectors clojure.lang.PersistentVector}
        [v idx & sxdi] (reverse idxsv)
        idxs (reverse sxdi)
        v (if-let [h (hints hint)] (list h v) v)
        nested-array (if (seq idxs)
                       `(deep-aget ~'objects ~array ~@idxs)
                        array)
        a-sym (with-meta (gensym "a") {:tag hint})]
      `(let [~a-sym ~nested-array]
         (aset ~a-sym ~idx ~v))))


;;; turn on box warnings
(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- init-matrix
  "We use Java arrays. Which means we are in mutation land, which we
  trade for speed. Also, the arrays are strictly contained within a
  single alignment computation (each alignment allocates their own
  array)."
  [rows cols endsgap]
  (let[scmat (make-array (class []) rows cols)
       endsgap (long endsgap)]
    (dotimes [i (long rows)]
      (dotimes [j (long cols)]
        (deep-aset objects scmat i j
              (cond (= i j 0) [:- 0 [ 0 0]]
                    (= i 0) [:l (* j endsgap) [0 j]]
                    (= j 0) [:u (* i endsgap) [i 0]]
                    :else []))))
    scmat))

(defn- score
  "Compute and return the best score and its direction :l for 'from
  left' :u for 'from up' and :d for 'from diagonal'. There is a lot of
  ugly stuff here to keep cpu cycles low. For example, we don't need
  to use an explicit substitution matrix and can use match and
  m(is)match scores directly (a 10x speed up) and we don't explicitly
  sort the resulting scores instead opting for a nested (if .. then
  else) approach. And a lot of primitive typing to ensure unboxed
  arithmatic. Also use of str/get removes chatAt reflection!! (or could
  have added ^String hints to s1 and s2, but ugly)"
  [i j c1 c2 scmat kind, gap match mmatch submat]
  ;;(println :I i :J j (map vec (vec scmat)))
  (let [i (long i)
        j (long j)
        gap (long gap)
        match (long match)
        mmatch (long mmatch)
        ;;c1 (str/get s1 i)
        ;;c2 (str/get s2 j)
        submat? submat
        lsc (+ (long ((deep-aget objects scmat i (dec j)) 1)) gap)
        usc (+ (long ((deep-aget objects scmat (dec i) j) 1)) gap)
        dsc (+ (long ((deep-aget objects scmat (dec i) (dec j)) 1))
               (if submat?
                 (long (submat [c1 c2]))
                 (if (= c1 c2) match mmatch)))
        score (if (> lsc usc)
                (if (> lsc dsc)
                  [:l lsc]
                  [:d dsc])
                (if (> usc dsc)
                  [:u usc]
                  [:d dsc]))
        ;; Check if local for possible start location of [i j]
        score (if (and (= kind :local) (> 0 (long (score 1)))) [:s 0] score)]
    (conj score [i j])))

(defmacro do-col [scmat col e & body]
  `(dotimes [r# (alength ~scmat)]
     (let [~e (deep-aget ~'objects ~scmat r# ~col)]
       ~@body)))

(defmacro do-row [scmat row e & body]
  `(let [scr# ^"[Lclojure.lang.PersistentVector;" (aget ~scmat 0)]
     (dotimes [c# (alength scr#)]
       (let [~e (deep-aget ~'objects ~scmat ~row c#)]
         ~@body))))

(defn- find-start
  "Find the starting cell in the score matrix for trace back. This is
  also the cell with the 'best' score for each KIND of alignment:

  :global - Standard NW global alignment
  :ends-gap-free - global alignment without gaps on ends (prefix/suffix aln)
  :local - Standard SW local alignment

  This also contains a lot of ugly low level matrix access code to
  avoid consing and vector overhead. Made somewhat better by do-col
  and do-row macros.
  "
  [^"[[Lclojure.lang.PersistentVector;" scmat, rows cols kind]
  (case kind
    :global (aget scmat rows cols)
    :ends-gap-free
    (let [max (volatile! (aget scmat rows cols))]
      (do-row scmat rows v
              (when (> (long (v 1)) (long (@max 1))) (vswap! max (fn[_] v))))
      (do-col scmat cols v
              (when (> (long (v 1)) (long (@max 1))) (vswap! max (fn[_] v))))
      @max)
    :local
    (let [max (volatile! (aget scmat rows cols))]
      (dotimes [r rows]
        (do-row scmat r v
                (when (> (long (v 1)) (long (@max 1))) (vswap! max (fn[_] v)))))
      @max)))

(defn- trace-back [scmat rows cols kind]
  (let [start (find-start scmat rows cols kind)]
    #_(clojure.pprint/pprint (map vec (vec scmat)))
    (loop [P (list start)
           [r c] (start 2)]
      #_(println r c P)
      (cond
        (and (= kind :global) (= r c 0)) (rest P)
        (and (not= kind :global) (= ((aget scmat r c) 1) 0)) (rest P)
        :else
        (let [r (long r)
              c (long c)
              cell (aget scmat r c)
              dir (cell 0)
              step (case dir
                     :d (deep-aget objects scmat (dec r) (dec c))
                     :u (deep-aget objects scmat (dec r) c)
                     :l (deep-aget objects scmat r (dec c))
                     (aum/raise
                      :dash-dir? "Bad direction"
                      :cell cell :r r :c c))]
          (recur (conj P step) (step 2)))))))

(defn decode-trace-back
  [s1 s2 tbk]
  (let [first-cell (first tbk)
        last-cell (last tbk)
        scr (second last-cell)
        start-idx (->> first-cell last (mapv dec))
        last-idx  (->> last-cell last (mapv dec))
        char-pair (fn[[dir scr [i j]]]
                    (case dir
                      :d [(str/get s1 i) (str/get s2 j)]
                      :l ["-" (str/get s2 j)]
                      :u [(str/get s1 i) "-"]))
        pairs (map char-pair tbk)]
    [[scr start-idx last-idx]
     [(apply str (map first pairs))
      (apply str (map second pairs))]]))

(defn align
  "Align s1 and s2 (strings) in the manner defined by kind:

  :global - Standard NW global alignment
  :ends-gap-free - global aln without gaps on ends (prefix/suffix aln)
  :local - Standard SW alignment

  submat is a substitution map (matrix) giving a score for a character
  pair. This namespace provides nt4submat which encodes a 'standard'
  ATGCxATGC substitution matrix.

  match is a score for the case when characters match (they are =)
  mmatch is a penalty score for the case when characters do not match

  submat and match/mmmatch are mutually exclusive

  gap is a penalty score for indels (characters match a gap)

  keep-scmat is a boolean for whether to return the scoring
  matrix (true) or not (false)

  Returns the traceback path of best alignment as a vector of
  triples: [dir score [i j]], where dir is the direction to get this
  point of path (:d -> diagonal, :u -> up, :l -> left), score is the
  score at this point, and [i j] is the index pair for s1 and s2 of
  this point (i for s1, j for s2) 0-based.

  If decode is true (default), tb is returned as the pair of aligned
  strings (with gaps if needed)

  If keep-scmat is true, returns a vector [tb scmat] otherwise returns
  just tb, tb the traceback.
  "
  [s1 s2 & {:keys [kind submat gap match mmatch decode keep-scmat]
            :or {kind :global gap -2 decode true keep-scmat false}}]
  (let [rows (count s1)
        cols (count s2)
        s1 (str "-" s1)
        s2 (str "-" s2)
        endsgap (if (= kind :global) gap 0)
        scmat (init-matrix (inc rows) (inc cols) endsgap)]
    (dotimes [c (long cols)]
      (dotimes [r (long rows)]
        (let [c (inc c)
              r (inc r)
              ch1 (str/get s1 r)
              ch2 (str/get s2 c)]
          (deep-aset objects scmat r c
                (score r c ch1 ch2 scmat kind gap match mmatch submat)))))
    (let [tb (trace-back scmat rows cols kind)
          tb (if decode (decode-trace-back s1 s2 tb) tb)]
      (if keep-scmat [tb scmat] tb))))


(defn- lcseq-score
  [i j s1 s2 scmat]
  (let [i (long i)
        j (long j)
        c1 (char (str/get s1 i))
        c2 (char (str/get s2 j))
        score (if (= c1 c2)
                [:d (+ 1 (long ((deep-aget objects scmat (dec i) (dec j)) 1)))]
                (let [usc (long ((deep-aget objects scmat (dec i) j) 1))
                      lsc (long ((deep-aget objects scmat i (dec j)) 1))]
                  (if (> lsc usc)
                    [:l lsc]
                    [:u usc])))]
    (conj score [i j])))

(defn lcsubseq
  [s1 s2]
  (let [rows (count s1)
        cols (count s2)
        s1 (str "-" s1)
        s2 (str "-" s2)
        scmat (init-matrix (inc rows) (inc cols) 0)]
    (dotimes [c (long cols)]
      (dotimes [r (long rows)]
        (let [c (inc c)
              r (inc r)]
          (deep-aset objects scmat r c
                (lcseq-score r c s1 s2 scmat)))))
    (trace-back scmat rows cols :global)))


(defn- lcstg-score
  [i j s1 s2 scmat]
  (let [i (long i)
        j (long j)
        c1 (char (str/get s1 i))
        c2 (char (str/get s2 j))
        score (if (= c1 c2)
                (+ 1 (long ((deep-aget objects scmat (dec i) (dec j)) 1)))
                0)]
    [:d score [i j]]))

(defn lcsubstg
  [s1 s2]
  (let [rows (count s1)
        cols (count s2)
        s1 (str "-" s1)
        s2 (str "-" s2)
        scmat (init-matrix (inc rows) (inc cols) 0)]
    (dotimes [c (long cols)]
      (dotimes [r (long rows)]
        (let [c (inc c)
              r (inc r)]
          (deep-aset objects scmat r c
                (lcstg-score r c s1 s2 scmat)))))
    (trace-back scmat rows cols :local)))


;;; Hirschberg's algorithm for linear space alignment
;;;
(defn- init-hirsch-matrix
  "We use Java arrays to mutate both cells and flip the rows back and
  forth. Arrays are strictly contained within a single alignment
  computation - each alignment (including subalignments of the divide
  and conquer) allocates their own array."
  [^long cols, ^long endsgap]
  (let [rows 2
        m (make-array (class []) rows cols)]
    (dotimes [i (long rows)]
      (dotimes [j (long cols)]
        (deep-aset objects m i j
                   (cond (= i j 0) [:- 0 [ 0 0]]
                         (= i 0) [:l (* j endsgap) [0 j]]
                         (= j 0) [:u (* i endsgap) [i 0]]
                         :else []))))
    m))

(defn- flip-rows
  "The matrix for Hirschberg's algo, is two rows (last two computed
  rows of full matrix). Only the second row is actively computed - the
  other is the original initial row (gap set) or the previous computed
  row. So, we only need two total rows and can (since we are in
  mutable Java land) reuse the previous computed row for the current
  to be computed row. So, we can merely interchange the rows. This
  simple function does that."
  [^"[[Lclojure.lang.PersistentVector;" hirsch-matrix]
  (let [m hirsch-matrix
        r0 (aget m 0)]
    (aset m 0 (aget m 1))
    (aset m 1 r0)
    m))

(defn- hirsch-score-row
  [s1 s2 & {:keys [kind submat gap match mmatch]}]
  (let [rows (count s1)
        cols (count s2)
        gap (long gap)
        s1 (str "-" s1)
        s2 (str "-" s2)
        ^"[[Lclojure.lang.PersistentVector;"
        rmat (init-hirsch-matrix (inc cols) gap)]
    (dotimes [r (long rows)]
      ;;(clojure.pprint/pprint rmat)
      (let [r (inc r)
            ch1 (str/get s1 r)]
        (dotimes [c (long cols)]
          (let [c (inc c)
                ch2 (str/get s2 c)]
            (deep-aset objects rmat 1 c
                  (score 1 c ch1 ch2 rmat :global gap match mmatch submat))))
        (flip-rows rmat)
        (deep-aset objects rmat 1 0 [:u (* (inc r) gap) [(inc r) 0]])))
    (aget rmat 0)))

(defn- hirsch-comb
  [scr l r]
  ;;(println :L l :R r)
  (let [l (if (integer? (first @l)) (second @l) @l)
        r (if (integer? (first @r)) (second @r) @r)]
    [scr (coll/concatv l r)]))

(defn hirsch-align-recur
  [sx sy & {:keys [kind submat gap match mmatch]}]
  #_(println :SX sx :SY sy)
  (cond
    (= (count sx) 0) [[(str/repeat (count sy) "-") sy]]
    (= (count sy) 0) [[sx (str/repeat (count sx) "-")]]
    (= (count sx) (count sy) 1) [[sx sy]]
    (or (= (count sx) 1) (= (count sy) 1))
    [(second (align sx sy :match match :mmatch mmatch :gap gap))]

    :else
    (let [mid (/ (count sx) 2)
          sxb (str/substring sx 0 mid)
          sxe (str/substring sx mid)
          S (future
              (hirsch-score-row sxb sy :match match :mmatch mmatch :gap gap))
          R (future
              (hirsch-score-row (str/reverse sxe) (str/reverse sy)
                                :match match :mmatch mmatch :gap gap))
          scl (->> @S vec (map (fn [[_ sc ij]] [sc ij])) (map first))
          scr (->> @R vec (map (fn [[_ sc ij]] [sc ij])) (map first))
          [ymid scr] (->> [scl (reverse scr)] (apply map +)
                          (keep-indexed #(vector %1 %2))
                          (apply max-key second))
          syb (str/substring sy 0 ymid)
          sye (str/substring sy ymid)]
      (hirsch-comb
       scr
       (future
         (hirsch-align-recur sxb syb :match match :mmatch mmatch :gap gap))
       (future
         (hirsch-align-recur sxe sye :match match :mmatch mmatch :gap gap))))))

#_(defn- hirsch-comb
  [scr l r]
  (println :SCR scr :L l :R r)
  (let [l (if (integer? (first l)) (second l) l)
        r (if (integer? (first r)) (second r) r)]
    [scr (coll/concatv l r)]))

#_(defn hirsch-align-recur
  [sx sy & {:keys [kind submat gap match mmatch]}]
  ;;(println :SX sx :SY sy)
  (cond
    (= (count sx) 0) [[(str/repeat (count sy) "-") sy]]
    (= (count sy) 0) [[sx (str/repeat (count sx) "-")]]
    (= (count sx) (count sy) 1) [[sx sy]]
    (or (= (count sx) 1) (= (count sy) 1))
    [(second (align sx sy :match match :mmatch mmatch :gap gap))]

    :else
    (let [mid (/ (count sx) 2)
          sxb (str/substring sx 0 mid)
          sxe (str/substring sx mid)
          S (hirsch-score-row sxb sy :match match :mmatch mmatch :gap gap)
          R (hirsch-score-row (str/reverse sxe) (str/reverse sy)
                              :match match :mmatch mmatch :gap gap)
          scl (->> S vec (map (fn [[_ sc ij]] [sc ij])) (map first))
          scr (->> R vec (map (fn [[_ sc ij]] [sc ij])) (map first))
          ;;_ (println :SCL scl :SCR scr)
          [ymid scr] (->> [scl (reverse scr)] (apply map +)
                          (keep-indexed #(vector %1 %2))
                          (apply max-key second))
          ;;_ (println :YMID ymid :SCORE scr)
          syb (str/substring sy 0 ymid)
          sye (str/substring sy ymid)]
      (hirsch-comb
       scr
       (hirsch-align-recur sxb syb :match match :mmatch mmatch :gap gap)
       (hirsch-align-recur sxe sye :match match :mmatch mmatch :gap gap)))))

(defn hirsch-align
  [sx sy & {:keys [kind submat gap match mmatch with-score? score-only?]
            :or {kind :global gap -2 with-score? true}}]
  (if score-only?
    (-> (hirsch-score-row sx sy :match match :mmatch mmatch :gap gap)
        last second)
    (let [[sc aln] (hirsch-align-recur
                    sx sy :match match :mmatch mmatch :gap gap)
          aln (reduce (fn[[X Y] [x y]] [(str X x) (str Y y)]) aln)]
      (if with-score? [sc aln] aln))))



;; Turn off box warnings - NOTE: this will be reset at end of file
;; anyway, but to be clear (documented) explicitly turn off
(set! *unchecked-math* false)
(set! *warn-on-reflection* false)


(comment

  (r/fold
   100
   (fn([] ["" ""])
     ([[L R] [ls rs]] [(apply str L ls) (apply str R rs)]))
   (fn[[L R] [l r]] [(str L l) (str R r)])
   (hirsch-align-recur
    sx sy :match match :mmatch mmatch :gap gap))

  (r/fold
   25
   (fn([] ["" ""])
     ([[L R] [ls rs]] [(apply str L ls) (apply str R rs)]))
   (fn[[L R] ls]
     [(apply str L (map first ls)) (apply str R (map second ls))])
   (coll/partitionv-all
    100
    (hirsch-align-recur
     sx sy :match match :mmatch mmatch :gap gap)))

  (hirsch-align "AGTACGCA" "TATGC" :match 2 :mmatch -1 :gap -2)


  (def tstsq "CTTCTGTCTATACCATTCTTCAAACCCTATAGTACATAAGTTTCCTTCTATAAATAGCTCGAAAAAAATTTATTTTTATGCCTGTCTCTTATACACATCTCCGAGCCCACGAGACAGGCAGAAATCTCGTATGCCGTCTTCTGCTTGAAA")

  (let [s1 tstsq
        s2 "CTGTCTCTTATACACATCT" ; transposon
        ;;s2 "ATCTCGTATGCCGTCTTCTGCTTG" ; p7 primer
        rows (count s1)
        cols (count s2)]
    (align s1 s2 :kind :ends-gap-free :match 2 :mmatch -1 :gap -2))

  (let [s1 "aaabca"
        s2 "abc"]
    (align s1 s2 :kind :ends-gap-free :gap -3 :match 5 :mmatch -2))
  )
