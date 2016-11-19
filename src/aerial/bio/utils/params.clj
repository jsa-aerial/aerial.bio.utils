;;--------------------------------------------------------------------------;;
;;                                                                          ;;
;;                   B I O . U T I L S . C O N F I G                        ;;
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

(ns aerial.bio.utils.params
  (:require
   [clojure.tools.reader.edn :as edn]
   [aerial.fs :as fs]))


;;; ------------------------------------------------------------------------;;;
;;;

(def params (atom {}))

(defn param-ks
  "Returns the set of keys indexing the params db. Only the first
  level"
  []
  (keys @params))

(defn param-set [] @params)

(defn get-params
  "Returns the value(s) indexed by ks. ks is either a single keyword,
  an enumeration of keywords, or a vector of keywords. The single
  keyword returns the top level params db value indexed by the
  keyword. An enumerated listing of keywords returns all top level
  params db values indexed by each keyword. A vector of keywords
  returns the multi-level value associated with the keyword path in
  the params db (like get-in)."
  [& ks]
  (if (= 1 (count ks))
    (let [k (first ks)]
      (if (vector? k)
        (get-in @params k)
        (@params k)))
    (mapv @params ks)))


;;; We need a fully general purpose way of setting default genome
;;; fasta directories as several functions use this in large call
;;; trees.  By making the default a function, we can have the
;;; 'directory' be multiple directories and different for different
;;; entries.

;;; Can't use delay here as these are intended to be dynamically
;;; bindable as well as set at configuration time.
;;;
(def ^{:dynamic true
       :doc "Default genome fasta location resolver"}
  default-genome-fasta-dir
  (get-params [:genomes :default]))

(def ^{:dynamic true
       :doc "Default blast database"}
  default-blast-db
  (fs/join (get-in @params [:blast :base])
           (get-in @params
                   [:blast (get-in @params [:blast :default])])))


(defn set-configuration
  "Sets up the params db and sets the default values for genome fasta
  dir and blast database. The basic shape of the config-map requires
  a :genomes submap and a :blast submap. Each of these requires
  a :default field and the value of the default field should be a
  keyword indexing the containing map (genomes/blast) and its value
  will be this indexed value.

  The config-map is presumed to have been obtained by reading a
  configuration map in a config file. For example, as by:

  (->> config-file-spec slurp clojure.tools.reader.edn/read-string)"
  [config-map]
  (let [genomes (config-map :genomes)
        base (genomes :base)
        default (genomes :default)
        genomes (->> (dissoc genomes :base :default)
                     (map (fn[[k v]] [k (fn[& _] (fs/join base v))]))
                     (into {})
                     (#(assoc % :base base :default (% default))))]
    (swap! params (fn[_] (assoc config-map :genomes genomes)))
    (alter-var-root
     #'default-genome-fasta-dir (fn[_] (get-params [:genomes :default])))
    (alter-var-root
     #'default-blast-db
     (fn[_]
       (fs/join (get-params [:blast :base])
                (get-params [:blast (get-params [:blast :default])]))))))


(defmacro with-genome-db
  "Binds DBREF as the current genome fasta directory resolver.  DBREF
   is a genome db fasta directory resolver, which is a function of 1
   argument returning the director for the genome name passed in.
   Executes BODY (a group of forms) within context so bound -
   including across threads
  "
  [dbref & body]
  `(with-bindings*
     {#'default-genome-fasta-dir ~dbref}
     (fn[] ~@body)))

(defmacro with-blast-db
  "Binds DBREF as the current blast database.  DBREF is a blast
  database specification as expected by NCBI blast cmds.  Executes
  BODY (a group of forms) within context so bound - including across
  threads
  "
  [dbref & body]
  `(with-bindings*
     {#'default-blast-db ~dbref}
     (fn[] ~@body)))
