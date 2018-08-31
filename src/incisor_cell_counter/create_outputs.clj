(ns incisor-cell-counter.create-outputs
  (:gen-class)
  (:import (mosaic.region_competition.RC AlgorithmRC)
           (mosaic.region_competition.DRS SobelVolume)
           (ij.process StackStatistics)
           (mosaic.core.imageUtils.images IntensityImage)
           (mosaic.plugins Region_Competition)
           (net.imglib2.view Views)
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

(def basename "/Users/kharrington/Data/Sharir_Amnon/Results005/")
(def counts (ij/open-img (str basename "24h-242Rc_cellCounts.tif")))
(def signals (ij/open-img (str basename "24h-242Rc_cellSignals.tif")))
(def membrane (ij/open-img (str basename "24h-242Rc_membrane.tif")))
(def sig (ij/open-img (str basename "24h-242Rc_signal.tif")))

(def min-count 10)
(def max-count (* 7 7 10))

(def tmp-001 (fun.imagej.ops.threshold/apply counts (net.imglib2.type.numeric.real.FloatType. (float max-count))))
(def tmp-002 (fun.imagej.ops.threshold/apply (fun.imagej.ops.math/multiply counts (net.imglib2.type.numeric.real.FloatType. (float -1)))
                                             (net.imglib2.type.numeric.real.FloatType. (float (- min-count)))))

(def size-filter-mask
  (fun.imagej.ops.convert/float32
  (fun.imagej.ops.math/multiply
    (let [maxthresh (fun.imagej.ops.threshold/apply counts (net.imglib2.type.numeric.real.FloatType. (float max-count)))]
      (fun.imagej.ops.image/invert (fun.imagej.ops.create/img maxthresh) maxthresh))
    (let [minthresh (fun.imagej.ops.threshold/apply (fun.imagej.ops.math/multiply counts (net.imglib2.type.numeric.real.FloatType. (float -1)))
                                                    (net.imglib2.type.numeric.real.FloatType. (float (- min-count))))]
      (fun.imagej.ops.image/invert (fun.imagej.ops.create/img minthresh) minthresh)))))
;(def size-filter-mask (fun.imagej.ops.convert/float32 size-filter-mask))

(def output (concat [sig membrane]
                    (doall (map img/normalize (map #(fun.imagej.ops.math/multiply size-filter-mask %)
                                                   [counts signals])))))
#_(def output (doall (map img/normalize (map #(fun.imagej.ops.math/multiply size-filter-mask %)
                                             [membrane sig (fun.imagej.ops.math/divide signals counts)]))))
(def out output)

(ij/save-img (second output) (str basename "24h-242Rc_counts.tif"))
(ij/save-img (last output) (str basename "24h-242Rc_signals.tif"))

(ij/show size-filter-mask (str "size filter mask [" min-count " to " max-count "]"))
(let [composite (ij.CompositeImage. (convert/img->imp (xyz-to-xycz out)))]
  (ij1/show-imp composite)
  (ij1/save-imp-as-tiff composite
                        (str basename "/24h-242Rc_composite.tif")))

;; up to here

  #_(ij1/show-imp (ij1/zconcat-imps
                  (map (fn [r g b]
                         (ij1/convert-stack-to-rgb (ij1/zconcat-imps (map convert/img->imp (map fun.imagej.ops.convert/uint8
                                                                                                (map #(fun.imagej.ops.math/multiply % (net.imglib2.type.numeric.real.FloatType. 255.0))
                                                                                                     [r g b]))))))
                       (img/dimension-split (first out) 2)
                       (img/dimension-split (second out) 2)
                       (img/dimension-split (last out) 2))))



  #_(ij1/show-imp (.mergeHyperstacks (ij.plugin.RGBStackMerge.)
                                   (into-array ij.ImagePlus (take 3 (map convert/img->imp (map fun.imagej.ops.convert/uint16 out)))) true))

  #_(ij1/show-imp (ij1/imps-to-rgb (map convert/img->imp out))); map dimension-split over channel-imgs (each image needs dimension-split, convert/img->imp, imps-to-rgb (per-slice), zconcat-imps)

  #_(ij1/show-imp (ij1/convert-to-RGB (convert/img->imp output)))
