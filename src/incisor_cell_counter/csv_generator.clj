(ns incisor-cell-counter.csv-generator
  (:gen-class)
  (:import (mosaic.region_competition.RC AlgorithmRC)
           (mosaic.region_competition.DRS SobelVolume)
           (ij.process StackStatistics)
           (mosaic.core.imageUtils.images IntensityImage)
           (mosaic.plugins Region_Competition)
           (net.imglib2.view Views)
           (net.imglib2 RandomAccessible RandomAccess)
           (loci.plugins.in ImporterOptions))
  (:require [fun.imagej.img :as img]
            [fun.imagej.imp :as ij1]
            [fun.imagej.conversion :as convert]
            [fun.imagej.img.cursor :as cursor]
            [fun.imagej.img.shape :as shape]
            [fun.imagej.img.type :as tpe]
            [fun.imagej.core :as ij]
    ;[fun.imagej.ops :as ops]
            [fun.imagej.img.utils :as img-utils]
            [fun.imagej.mesh :as msh]
            [fun.imagej.segmentation.fast-segmentation :as fseg]
            [clojure.java.io :as io]
            [clojure.string :as string]))

(def display false)

(when display
  (defonce ijui (ij/show-ui)))

(defn xyz-to-xycz
  "Convert a seq of xyz imgs into a xycz stack"
  [imgs]
  (let [dimensions (img/dimensions (first imgs))
        stack (fun.imagej.ops.create/img (long-array [(first dimensions) (second dimensions) (count imgs) (last dimensions)]))]
    (dotimes [z (last dimensions)]
      (dotimes [c (count imgs)]
        (img/map-img! cursor/copy-real
                      (img/hyperslice (img/hyperslice stack 3 z) 2 c); slice into C and Z
                      (img/hyperslice (nth imgs c) 2 z))))
    stack))

(defn xycz-to-xyz
  "Convert a seq of xyz imgs into a xycz stack"
  [stack]
  (let [dimensions (seq (img/dimensions stack))
        imgs (doall (for [k (range (nth dimensions 2))]
                      (fun.imagej.ops.create/img (long-array [(first dimensions) (second dimensions) (last dimensions)]))))]
    (dotimes [z (last dimensions)]
      (dotimes [c (count imgs)]
        (img/map-img! cursor/copy-real
                      (img/hyperslice (nth imgs c) 2 z)
                      (img/hyperslice (img/hyperslice stack 3 z) 2 c))))
    imgs))

(defn threshold-signal
  "Return a signal thresholded by 250% moments threshold"
  [target]
  (let [input (img/copy target)
        hist (fun.imagej.ops.image/histogram input)
        thresh-val (* 2.5 (.get (fun.imagej.ops.threshold/moments hist)))]
    (fun.imagej.ops.convert/uint8 (fun.imagej.ops.threshold/apply input (net.imglib2.type.numeric.real.FloatType. thresh-val)))))

(set! *warn-on-reflection* true)

(defn -main
  [& args]
  (println "Starting")
  (let [directory (or (first args)
                      "/Users/kharrington/Data/Sharir_Amnon/Results009/")
        just-name (or (second args)
                      "45m-4001c")
        basename (str directory just-name)]
  (def cell-labels (ij/open-img (str basename "_cellLabels.tif")))
  (def sig (ij/open-img (str basename "_signal.tif")))
  (def binary-signal (threshold-signal sig))
  (when display
    (ij/show binary-signal))
  (println "Starting counts")
  (def voxel-res {:w 1 :h 1 :d 2})

  (def voxel-counts (atom {}))
  (def signal-counts (atom {}))
  (def binary-signal-counts (atom {}))
  (def cell-xs (atom {}))
  (def cell-ys (atom {}))
  (def cell-zs (atom {}))

  (def crop-id (atom {}))
  (def min-xs (atom {}))
  (def min-ys (atom {}))
  (def min-zs (atom {}))

  (def max-xs (atom {}))
  (def max-ys (atom {}))
  (def max-zs (atom {}))

  (let [img-dims (seq (img/dimensions cell-labels))
        crop-dims [100 100 (last img-dims)]]
    ; Calculating stats
    ; Signal stats
    ; Size and position stats
    (let [pos ^longs (long-array 3)]
      (img/map-img! (fn [label binary-signal signal]
                      (.localize ^net.imglib2.Cursor label pos)
                      (when (and (zero? (aget pos 0))
                                 (zero? (aget pos 1)))
                        (println (seq pos)))
                      (let [start-x (int (* (first crop-dims) (java.lang.Math/floor (/ (aget pos 0) (first crop-dims)))))
                            start-y (int (* (second crop-dims) (java.lang.Math/floor (/ (aget pos 1) (second crop-dims)))))
                            lid (str start-x "_" start-y "_" (.get ^net.imglib2.type.numeric.real.FloatType (.get  ^net.imglib2.Cursor label)))]
                        ;(println lid)
                        (reset! signal-counts
                                (assoc @signal-counts
                                  lid (if (get @signal-counts lid)
                                        (+ (get @signal-counts lid) (.get ^net.imglib2.type.numeric.real.FloatType (.get ^net.imglib2.Cursor signal)))
                                        (.get ^net.imglib2.type.numeric.real.FloatType (.get ^net.imglib2.Cursor signal)))))
                        (reset! binary-signal-counts
                                (assoc @binary-signal-counts
                                  lid (if (get @binary-signal-counts lid)
                                        (+ (get @binary-signal-counts lid) (.get ^net.imglib2.type.numeric.integer.UnsignedByteType (.get  ^net.imglib2.Cursor binary-signal)))
                                        (.get ^net.imglib2.type.numeric.integer.UnsignedByteType (.get  ^net.imglib2.Cursor binary-signal)))))
                        (reset! voxel-counts
                                (assoc @voxel-counts
                                  lid (if (get @voxel-counts lid)
                                        (inc (get @voxel-counts lid))
                                        1)))
                        (reset! crop-id
                                (assoc @voxel-counts
                                  lid (if (get @voxel-counts lid)
                                        (get @voxel-counts lid)
                                        (.get ^net.imglib2.type.numeric.real.FloatType (.get  ^net.imglib2.Cursor label)))))
                        (reset! cell-xs
                                (assoc @cell-xs
                                  lid (if (get @cell-xs lid)
                                        (+ (aget pos 0) (get @cell-xs lid))
                                        (aget pos 0))))
                        (reset! cell-ys
                                (assoc @cell-ys
                                  lid (if (get @cell-ys lid)
                                        (+ (aget pos 1) (get @cell-ys lid))
                                        (aget pos 1))))
                        (reset! cell-zs
                                (assoc @cell-zs
                                  lid (if (get @cell-zs lid)
                                        (+ (aget pos 2) (get @cell-zs lid))
                                        (aget pos 2))))
                        (reset! min-xs
                                (assoc @min-xs
                                  lid (if (get @min-xs lid)
                                        (min (aget pos 0) (get @min-xs lid))
                                        (aget pos 0))))
                        (reset! min-ys
                                (assoc @min-ys
                                  lid (if (get @min-ys lid)
                                        (min (aget pos 1) (get @min-ys lid))
                                        (aget pos 1))))
                        (reset! min-zs
                                (assoc @min-zs
                                  lid (if (get @min-zs lid)
                                        (min (aget pos 2) (get @min-zs lid))
                                        (aget pos 2))))
                        (reset! max-xs
                                (assoc @max-xs
                                  lid (if (get @max-xs lid)
                                        (max (aget pos 0) (get @max-xs lid))
                                        (aget pos 0))))
                        (reset! max-ys
                                (assoc @max-ys
                                  lid (if (get @max-ys lid)
                                        (max (aget pos 1) (get @max-ys lid))
                                        (aget pos 1))))
                        (reset! max-zs
                                (assoc @max-zs
                                  lid (if (get @max-zs lid)
                                        (max (aget pos 2) (get @max-zs lid))
                                        (aget pos 2))))
                        ))
                    cell-labels binary-signal sig))
    (println "Writing stats to file")
    (let [ks (keys @voxel-counts)
          max-val (dec (java.lang.Math/pow 2 31))]
      (spit (str basename "_cellStats.csv")
            (with-out-str
              (println (string/join "," ["cellID" "x" "y" "z" "voxels" "signalSum" "min_x" "min_y" "min_z" "max_x" "max_y" "max_z" "binarySignal"]))
              (doseq [k ks]
                (when (and (< (get @crop-id k) max-val)
                           (> (get @voxel-counts k) 1)
                  (println (string/join "," [k
                                             (float (/ (get @cell-xs k) (get @voxel-counts k)))
                                             (float (/ (get @cell-ys k) (get @voxel-counts k)))
                                             (float (/ (get @cell-zs k) (get @voxel-counts k)))
                                             (get @voxel-counts k) (or (get @signal-counts k) 0)
                                             (get @min-xs k) (get @min-ys k) (get @min-zs k)
                                             (get @max-xs k) (get @max-ys k) (get @max-zs k)
                                             (get @binary-signal-counts k)])))))))))

  ;(def basename "/Users/kharrington/Data/Sharir_Amnon/Results007/")
  ;(def csv-name "_cellStats.csv")
  ;(def csv-name "24h-242Rc_cellStats_shape_filtered_firstpass.csv")

  (let [min-size 30
        max-size 300 ;5000
        min-width 2
        min-height 2
        min-depth 2]
    (defn segment-predicate
      "expecting a segment where leading element is id"
      [attrs]
      (let [width (- (nth attrs 9)
                     (nth attrs 6))
            height (- (nth attrs 10)
                      (nth attrs 7))
            depth (- (nth attrs 11)
                     (nth attrs 8))]
        ;(println width height depth attrs)
        (and (> (nth attrs 4) min-size)
             (< (nth attrs 4) max-size)
             (> width min-width)
             (> height min-height)
             (> depth min-depth)))))
  #_(take-while #(not (segment-predicate %))
                segments)
  ;(map segment-predicate (take 100 segments))

  (def filtered-segments
    (let [;lines (string/split-lines (slurp "/Users/kharrington/Google Drive/LeoMouseTeeth/withBinary/24h-242Rc_cellStats.csv"))
          ;lines (string/split-lines (slurp "/Users/kharrington/Google Drive/LeoMouseTeeth/24h-242Rc_cellStats.csv"))
          lines (string/split-lines (slurp (str basename "_cellStats.csv"))) ;/Users/kharrington/Google\ Drive/LeoMouseTeeth/withBinary/24h-242Rc_cellStats.csv
          header (string/split (first lines) #",")
          body (rest lines)
          segments (doall
                     (for [line body]
                       (let [parts (string/split line #",")
                             vs (map read-string (rest parts))
                             crop-parts (map read-string (string/split (first parts) #"_"))]
                         ;(println crop-parts (take 3 vs))
                         (concat [(first parts)
                                  (first vs)
                                  (second vs)
                                  #_(+ (first crop-parts) (first vs))
                                  #_(+ (second crop-parts) (second vs))]
                                 (drop 2 vs)))))
          filtered-segments (filter #_#(and (> (nth % 4) min-size)
                                            (< (nth % 4) max-size))
                              segment-predicate
                              segments)
          ]
      (def header header)
      (def segments segments)
      (zipmap (map first filtered-segments)
              (map rest filtered-segments))))
  (println :num-filtered-segments (count (keys filtered-segments)))
  (println (sort-by first (frequencies (map #(nth % 3) (vals filtered-segments)))))

  (def img-size (seq (img/dimensions cell-labels)))
  ;(def signals (threshold-signal (ij/open-img (str basename "24h-242Rc_signal.tif"))))
  ;(def signals (ij/open-img (str basename "24h-242Rc_signal.tif")))
  ;(def original-labels (ij/open-img (str basename "50crop/24h-242Rc_cellLabels.tif")))
  ;(def original-labels (ij/open-img (str basename "24h-242Rc_originalCellLabels.tif")))
  (def cell-labels (fun.imagej.ops.convert/uint16 (fun.imagej.ops.create/img (long-array img-size))))
  (def cell-signals (fun.imagej.ops.convert/float32 (fun.imagej.ops.create/img (long-array img-size))))
  (def avg-signals (fun.imagej.ops.convert/float32 (fun.imagej.ops.create/img (long-array img-size))))
  (def sphere-signals (fun.imagej.ops.convert/float32 (fun.imagej.ops.create/img (long-array img-size))))
  (def sphere-labels (fun.imagej.ops.convert/uint16 (fun.imagej.ops.create/img (long-array img-size))))

  ;; Previous masking
  #_(let [mem (ij/open-img (str basename "_membrane.tif"))]
      (def mask (let [full-mask (fun.imagej.ops.threshold/triangle ;previously intermodes
                                  (img/concat-imgs
                                    (doall (for [k (range (img/get-size-dimension mem 2))]
                                             (let [source (img/hyperslice mem 2 k)
                                                   output (fun.imagej.ops.create/img source)]
                                               (fun.imagej.ops.filter/gauss output source (double-array [2 2]))
                                               output
                                               #_(fun.imagej.ops.threshold/triangle output))))))]
                  (fun.imagej.ops.morphology/erode
                    full-mask (shape/sphere 10)))))



  (def mask (if (.exists (java.io.File. (str basename "_mergeMask.tif")))
              (fun.imagej.ops.threshold/apply (ij/open-img (str basename "_mergeMask.tif"))
                                            (net.imglib2.type.numeric.integer.UnsignedByteType. 1))
              (let [mem (ij/open-img (str basename "_membrane.tif"))]
                (let [full-mask (fun.imagej.ops.threshold/triangle ;previously intermodes
                                  (img/concat-imgs
                                    (doall (for [k (range (img/get-size-dimension mem 2))]
                                             (let [source (img/hyperslice mem 2 k)
                                                   output (fun.imagej.ops.create/img source)]
                                               (fun.imagej.ops.filter/gauss output source (double-array [2 2]))
                                               output
                                               #_(fun.imagej.ops.threshold/triangle output))))))]
                  (fun.imagej.ops.morphology/erode
                    full-mask (shape/sphere 10))))))
  (println "Mask loaded")

  ;(ij/save-img (fun.imagej.ops.convert/uint8 mask) "/Users/kharrington/Data/Sharir_Amnon/Results014/45m-0062d_initialMask.tif")

  ;(def include-mask (ij/open-img (str "/Users/kharrington/Data/Sharir_Amnon/Results014/45m-0062d_includeMask.tif")))
  ;(def exclude-mask (ij/open-img (str "/Users/kharrington/Data/Sharir_Amnon/Results014/45m-0062d_excludeMask.tif")))

  #_(def include-mask
    (ij/open-img (str basename "_includeMask.tif")))

  #_(def exclude-mask
    (ij/open-img (str basename "_excludeMask.tif")))


  ;(ij/show mask)

  (set! *warn-on-reflection* true)

  (println "processing")

  (def signal-spheres (atom {}))
  ; Draw spheres
  (let [pos ^longs (long-array 3)
        sphere-radius 3
        nbrhd-labels ^RandomAccess (.randomAccess (.neighborhoodsRandomAccessibleSafe (shape/sphere sphere-radius)
                                                                                      (Views/extendZero cell-labels)))
        nbrnd-signals ^RandomAccess (.randomAccess (.neighborhoodsRandomAccessibleSafe (shape/sphere sphere-radius)
                                                                                       (Views/extendZero cell-signals)))
        nbrnd-avg-signals ^RandomAccess (.randomAccess (.neighborhoodsRandomAccessibleSafe (shape/sphere sphere-radius)
                                                                                           (Views/extendZero avg-signals)))
        nbrnd-src-signal ^RandomAccess (.randomAccess (.neighborhoodsRandomAccessibleSafe (shape/sphere sphere-radius)
                                                                                          (Views/extendZero sig)))]
    (doseq [[seg-id attrs] filtered-segments]
      (dotimes [k 3]
        ;(println (nth attrs k))
        (aset pos k (long (nth attrs k))))
      (.setPosition nbrhd-labels pos)
      (.setPosition nbrnd-signals pos)
      (.setPosition nbrnd-avg-signals pos)
      (.setPosition nbrnd-src-signal pos)
      #_(println seg-id attrs (seq pos))
      (when (img/get-val mask pos)
        (let [nbr-label (.get nbrhd-labels)
              nbr-signal (.get nbrnd-signals)
              nbr-src-signal (.get nbrnd-src-signal)
              nbr-avg-signal (.get nbrnd-avg-signals)]
          (loop [label ^net.imglib2.Cursor (.cursor nbr-label)
                 signal ^net.imglib2.Cursor (.cursor nbr-signal)
                 avg-signal ^net.imglib2.Cursor (.cursor nbr-avg-signal)
                 src-signal ^net.imglib2.Cursor (.cursor nbr-src-signal)]
            (when (.hasNext label)
              (.fwd label) (.fwd signal) (.fwd avg-signal) (.fwd src-signal)
              (reset! signal-spheres
                      (assoc @signal-spheres
                        seg-id (if (get @signal-spheres seg-id)
                                 (+ (get @signal-spheres seg-id) (.get (.get src-signal)))
                                 (.get (.get src-signal)))))
              (.set (.get label) (hash seg-id))
              (.set (.get signal) (float (nth attrs 4)))
              (.set (.get avg-signal) (float (/ (nth attrs 4)
                                                (nth attrs 3))))
              (recur label signal avg-signal src-signal)))))))

  ; Drawing the sphere signal img
  (let [pos ^longs (long-array 3)
        cell-count (atom 1)
        sphere-radius 3
        nbrnd-sphere-labels  ^RandomAccess (.randomAccess (.neighborhoodsRandomAccessibleSafe (shape/sphere sphere-radius)
                                                                                              (Views/extendZero sphere-labels)))
        nbrnd-sphere-signals ^RandomAccess (.randomAccess (.neighborhoodsRandomAccessibleSafe (shape/sphere sphere-radius)
                                                                                              (Views/extendZero sphere-signals)))]
    (doseq [[seg-id attrs] filtered-segments]
      (dotimes [k 3]
        ;(println (nth attrs k))
        (aset pos k (long (nth attrs k))))
      (.setPosition nbrnd-sphere-signals pos)
      (.setPosition nbrnd-sphere-labels pos)
      (when (img/get-val mask pos)
        (let [nbr-sphere-signal (.get nbrnd-sphere-signals)
              nbr-sphere-label (.get nbrnd-sphere-labels)]
          (loop [sphere-signal ^net.imglib2.Cursor (.cursor nbr-sphere-signal)
                 sphere-label ^net.imglib2.Cursor (.cursor nbr-sphere-label)]
            (when (.hasNext sphere-signal)
              (.fwd sphere-signal)
              (.fwd sphere-label)
              (.set (.get sphere-signal) (float (get @signal-spheres seg-id)))
              (.set (.get sphere-label) (int @cell-count))
              (recur sphere-signal
                     sphere-label))))
        (swap! cell-count inc))))

  (println "done processing")

  #_(def signal-threshold
    (let [histogram (fun.imagej.ops.image/histogram sphere-signals)
          hist-val (fun.imagej.ops.threshold/renyiEntropy histogram)]
      (println "Signal cutoff for: " basename " is " hist-val (.get hist-val))
      hist-val))

  (def threshold-table
    {:45m 22500
     :24h 55000
     :48h 32000})

  (def signal-threshold
    (if (.contains basename "24h")
      (:24h threshold-table)
      (if (.contains basename "45m")
        (:45m threshold-table)
        (:48h threshold-table))))

  ; 45m 1 1 2
  ; 24h 1 1 2
  ; 48h 1 1 2

  (spit (str basename "_filteredSegments.csv")
        (with-out-str
          (println (string/join "," (concat header ["sampleSignal" "threshSignal"])))
          (doseq [k (keys filtered-segments)]
            (println (string/join "," (concat [k]
                                              (map * (take 3 (get filtered-segments k))
                                                   [(:w voxel-res) (:h voxel-res) (:d voxel-res)])
                                              (drop 3 (get filtered-segments k))
                                              [(or (get @signal-spheres k) 0)
                                               (if (and (get @signal-spheres k)
                                                        (> (get @signal-spheres k) signal-threshold))
                                                 1
                                                 0)]))))))

  (ij/save-img cell-labels (str basename "_filteredCellLabels.tif"))
  (ij/save-img cell-signals (str basename "_filteredCellSignals.tif"))
  (ij/save-img avg-signals (str basename "_filteredCellAvgSignals.tif"))
  (ij/save-img sphere-signals (str basename "_filteredCellSignalSamples.tif"))
  (ij/save-img (fun.imagej.ops.convert/uint8 (fun.imagej.ops.threshold/apply sphere-signals (net.imglib2.type.numeric.real.FloatType. (float signal-threshold))))
               (str basename "_filteredCellSignalSamples_binary.tif"))
  (ij/save-img sphere-labels (str basename "_filteredCellLabelSpheres.tif"))
  ))

;(-main)

;; Computing the mask from membrane
    ;    full-mask (fun.imagej.ops.threshold/triangle ;previously intermodes
    ;                (img/concat-imgs
    ;                  (doall (for [k (range (img/get-size-dimension mem 2))]
    ;                           (let [source (img/hyperslice mem 2 k)
    ;                                 output (fun.imagej.ops.create/img source)]
    ;                             (fun.imagej.ops.filter/gauss output source (double-array [2 2]))
    ;                             output
    ;                             #_(fun.imagej.ops.threshold/triangle output))))))]
    ;(fun.imagej.ops.morphology/erode
    ;  full-mask (shape/sphere 6))))

; For testing the mask
(do
  ;(def mem (ij/open-img "/Users/kharrington/Data/Sharir_Amnon/Results013/45m-0062d_membrane.tif"))
  ;(def mem (ij/open-img "/Users/kharrington/Data/Sharir_Amnon/Results014/48h-003d_membrane.tif"))
  ;(def mem (ij/open-img "/Users/kharrington/Data/Sharir_Amnon/Results015/48-003z_membrane.tif"))
  (def mem (ij/open-img "/Users/kharrington/Data/Sharir_Amnon/Results015/24h-242az_membrane.tif"))
  (def mask (let [full-mask (fun.imagej.ops.threshold/triangle ;previously intermodes
                                (img/concat-imgs
                                  (doall (for [k (range (img/get-size-dimension mem 2))]
                                           (let [source (img/hyperslice mem 2 k)
                                                 output (fun.imagej.ops.create/img source)]
                                             (fun.imagej.ops.filter/gauss output source (double-array [2 2]))
                                             output
                                             #_(fun.imagej.ops.threshold/triangle output))))))]
                (fun.imagej.ops.morphology/erode
                  full-mask (shape/sphere 10))))
  (ij/save-img (fun.imagej.ops.convert/uint8 mask) "/Users/kharrington/Data/Sharir_Amnon/Results015/24h-242az_initialMask.tif")
  )
