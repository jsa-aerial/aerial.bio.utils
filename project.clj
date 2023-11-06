(defproject aerial.bio.utils "2.2.1"
  :description "Various bioinformatic general utility tools"
  :url "https://github.com/jsa-aerial/aerial.bio.utils"
  :license {:name "The MIT License (MIT)"
            :url  "http://opensource.org/licenses/MIT"
            :distribution :repo}
  :dependencies
  [[org.clojure/clojure "1.10.3"]
   [org.clojure/data.csv "1.0.0"]
   [org.clojure/tools.reader "1.3.6"]

   [aerial.fs "1.1.6"]
   [aerial.utils "1.2.0"]

   ;; The following pulls in dnsjava, but it is required in utils so,
   ;; this should be pulled in by utils and dnsjava from that. But it
   ;; is not, and I as yet cannot figure out why.
   [net.apribase/clj-dns "0.1.0"]]

  :scm {:name "git"
        :url "https://github.com/jsa-aerial/aerial.bio.utils"}
  )
