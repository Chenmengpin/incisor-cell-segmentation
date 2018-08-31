(ns incisor-cell-counter.csv-to-image
  (:gen-class)
  (:import (mosaic.region_competition.RC AlgorithmRC)
           (mosaic.region_competition.DRS SobelVolume)
           (ij.process StackStatistics)
           (mosaic.core.imageUtils.images IntensityImage)
           (mosaic.plugins Region_Competition)
           (net.imglib2.view Views)
           (loci.plugins.in ImporterOptions)
           (net.imglib2 RandomAccessible RandomAccess))
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

(def display true)

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

(def basename "/Users/kharrington/Data/Sharir_Amnon/Results006/")
(def csv-name "24h-242Rc_cellStats_shape_filtered_firstpass.csv")
(def filtered-segments
  (let [lines (string/split-lines (slurp (str basename csv-name)))
        header (first lines)
        body (rest lines)
        segments (doall
                  (for [line body]
                    (let [vs (butlast ; drop extra filepath colm
                              (rest ; rest to drop the extra column
                                (map read-string (string/split line #","))))
                          crop-parts (map read-string (string/split (first vs) #"_"))]
                      (concat [(first vs)
                               (+ (second vs) (first crop-parts))
                               (+ (nth vs 2) (second crop-parts))]
                              (drop 3 vs)))))]
    (zipmap (map first segments)
            (map rest segments))))

(def img-size [1291 1192 175])
;(def original-labels (ij/open-img (str basename "50crop/24h-242Rc_cellLabels.tif")))
(def original-labels (ij/open-img (str basename "24h-242Rc_originalCellLabels.tif")))
(def cell-labels (fun.imagej.ops.convert/uint16 (fun.imagej.ops.create/img (long-array img-size))))
(def cell-signals (fun.imagej.ops.convert/float32 (fun.imagej.ops.create/img (long-array img-size))))

(set! *warn-on-reflection* true)

(println "processing")

(let [pos ^longs (long-array 3)
      nbrhd-labels ^RandomAccess (.randomAccess (.neighborhoodsRandomAccessibleSafe (shape/sphere 5)
                                                                                    (Views/extendZero cell-labels)))
      nbrnd-signals ^RandomAccess (.randomAccess (.neighborhoodsRandomAccessibleSafe (shape/sphere 5)
                                                                                     (Views/extendZero cell-signals)))]
  (doseq [[seg-id attrs] filtered-segments]
    (dotimes [k 3]
      ;(println (nth attrs k))
      (aset pos k (long (nth attrs k))))
    (.setPosition nbrhd-labels pos)
    (.setPosition nbrnd-signals pos)
    (let [nbr-label (.get nbrhd-labels)
          nbr-signal (.get nbrnd-signals)]
      (loop [label (.cursor nbr-label)
             signal (.cursor nbr-signal)]
        (when (.hasNext label)
          (.fwd label) (.fwd signal)
          (.set (.get label) (hash seg-id))
          (.set (.get signal) (float (nth attrs 4))))))))

; Fill using label img
#_(let [pos ^longs (long-array 3)]
  (img/map-img!
    (fn [^net.imglib2.Cursor original ^net.imglib2.Cursor label ^net.imglib2.Cursor signal]
      (.localize ^net.imglib2.Cursor original pos)
      (let [start-x (java.lang.Math/floor (aget pos 0))
            start-y (java.lang.Math/floor (aget pos 1))
            lid (str start-x "_" start-y "_" (float (.get ^net.imglib2.type.numeric.integer.UnsignedShortType (.get ^net.imglib2.Cursor original))))]
        ;(println lid (get filtered-segments lid))
        (when-let [props (get filtered-segments lid)]
          (.set ^net.imglib2.type.numeric.integer.UnsignedByteType (.get label) (int (.get ^net.imglib2.type.numeric.real.FloatType (.get original))))
          (.set ^net.imglib2.type.numeric.real.FloatType (.get signal) (float (nth props 5))))))
    original-labels cell-labels cell-signals))

(println "done processing")

(ij/save-img cell-labels (str basename "24h-242Rc_filteredCellLabels.tif"))
(ij/save-img cell-signals (str basename "24h-242Rc_filteredCellSignals.tif"))
