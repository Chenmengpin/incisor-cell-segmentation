(ns incisor-cell-counter.stitch-crops
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

(defn load-imaris-rl2
  [filename]
  (let [options (ImporterOptions.)
        _ (do (.clearSeries options)
              (.setId options filename)
              (.setSplitChannels options true)
              (.setSeriesOn options 1 true))
        im (loci.plugins.BF/openImagePlus options)]
    im
    ;im (open-bioformats "/Users/kharrington/Box Sync/Cervical loop kinetics -AS (Amnon Sharir)/Imaris images/48 hr/482-48h.ims")
    #_(ij1/show-imp (first im))))

(defn threshold-signal
  "Return a signal thresholded by 250% moments threshold"
  [target]
  (let [input (img/copy target)
        hist (fun.imagej.ops.image/histogram input)
        thresh-val (* 2.5 (.get (fun.imagej.ops.threshold/moments hist)))]
    (fun.imagej.ops.convert/uint8 (fun.imagej.ops.threshold/apply input (net.imglib2.type.numeric.real.FloatType. thresh-val)))))

(defn -main
  [& args]
  (println "Starting")
  (println "args" args)
  (let [directory (or (first args)
                      "/Users/kharrington/Data/Sharir_Amnon/Results007/")
        just-name (or (second args)
                      "24h-242Rc")
        basename (str directory just-name)

        [img-w img-h img-d crop-w crop-h] (if (> (count args) 6)
                                            (map #(read-string %) (take 5 (drop 2 args)))
                                            [100 200 175 50 50])
        _ (println [img-w img-h img-d crop-w crop-h])
        membrane (fun.imagej.ops.create/img (long-array [img-w img-h img-d]))
        signal (fun.imagej.ops.create/img (long-array [img-w img-h img-d]))
        cell-labels (fun.imagej.ops.create/img (long-array [img-w img-h img-d]))
        cell-signals (fun.imagej.ops.create/img (long-array [img-w img-h img-d]))
        cell-binary-signals (fun.imagej.ops.create/img (long-array [img-w img-h img-d]))
        cell-counts (fun.imagej.ops.create/img (long-array [img-w img-h img-d]))

        voxel-counts (atom {})
        signal-counts (atom {})
        binary-signal-counts (atom {})
        signal-val-lists (atom {})
        cell-xs (atom {})
        cell-ys (atom {})
        cell-zs (atom {})

        crop-id (atom {})
        min-xs (atom {})
        min-ys (atom {})
        min-zs (atom {})

        max-xs (atom {})
        max-ys (atom {})
        max-zs (atom {})

        img-dims [img-w img-h img-d]
        crop-dims [crop-w crop-h (last img-dims)]]
    (doall
      (for [start-x (range 0 (first img-dims) (first crop-dims))
            start-y (range 0 (second img-dims) (second crop-dims))]
        (let [start-point [start-x start-y 0]
              _ (println "Starting crop " start-point)
              offset (long-array start-point)
              stop-point (long-array (map #(dec (min (+ %1 %2)
                                                     %3))
                                          offset
                                          crop-dims
                                          img-dims))
              ; cropping
              membrane (net.imglib2.view.Views/offset (net.imglib2.view.Views/interval membrane offset stop-point)
                                                      (long-array start-point))
              signal (net.imglib2.view.Views/offset (net.imglib2.view.Views/interval signal offset stop-point)
                                                    (long-array start-point))
              cell-labels (net.imglib2.view.Views/offset (net.imglib2.view.Views/interval cell-labels offset stop-point)
                                                         (long-array start-point))
              cell-signals (net.imglib2.view.Views/offset (net.imglib2.view.Views/interval cell-signals offset stop-point)
                                                          (long-array start-point))
              cell-counts (net.imglib2.view.Views/offset (net.imglib2.view.Views/interval cell-counts offset stop-point)
                                                         (long-array start-point))
              ;inverted (normalize-in-mask inverted mask)
              crop-name (str "crop_sX_" start-x "_crop_sY_" start-y)

              _ (println (str basename "_" crop-name "_filtered.tif"))
              filtered-regions (ij/open-img (str basename "_" crop-name "_filtered.tif"))

              crop-imgs (xycz-to-xyz filtered-regions)]
          ; Copying images
          (img/map-img! cursor/copy cell-labels (first crop-imgs))
          (img/map-img! cursor/copy cell-signals (second crop-imgs))
          (img/map-img! cursor/copy cell-counts (nth crop-imgs 2))
          (img/map-img! cursor/copy membrane (nth crop-imgs 3))
          (img/map-img! cursor/copy signal (last crop-imgs))
          (println start-x start-y (img/get-val (first crop-imgs)
                                                (long-array [25 25 100])))
          ; Calculating stats
          ; Signal stats
          #_(img/map-img! (fn [label signal]
                          (let [lid (str start-x "_" start-y "_" (.get (.get label)))]
                            (swap! signal-val-lists
                                   (assoc @signal-val-lists
                                     lid (if (get @signal-val-lists lid)
                                           (conj (get @signal-val-lists lid) (.get (.get signal)))
                                           [(.get (.get signal))]))
                                   signal-val-lists)
                            (reset! signal-counts
                                    (assoc @signal-counts
                                      lid (if (get @signal-counts lid)
                                            (+ (get @signal-counts lid) (.get (.get signal)))
                                            (.get (.get signal)))))))
                        cell-labels signal)
          ; Size and position stats
          #_(let [pos ^longs (long-array 3)]
            (img/map-img! (fn [label]
                            (let [lid (str start-x "_" start-y "_" (.get (.get label)))]
                              (.localize label pos)
                              (reset! voxel-counts
                                      (assoc @voxel-counts
                                        lid (if (get @voxel-counts lid)
                                              (inc (get @voxel-counts lid))
                                              1)))
                              (reset! crop-id
                                      (assoc @voxel-counts
                                        lid (if (get @voxel-counts lid)
                                              (get @voxel-counts lid)
                                              (.get (.get label)))))
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
                          cell-labels)))))
    ;; Rerun on binary signal
    #_(let [binary-signal (threshold-signal signal)]
      (ij/save-img (fun.imagej.ops.convert/float32 binary-signal) (str basename "_binarySignal.tif"))
      (doall
        (for [start-x (range 0 (first img-dims) (first crop-dims))
              start-y (range 0 (second img-dims) (second crop-dims))]
          (let [start-point [start-x start-y 0]
                _ (println "Starting crop " start-point)
                offset (long-array start-point)
                stop-point (long-array (map #(dec (min (+ %1 %2)
                                                       %3))
                                            offset
                                            crop-dims
                                            img-dims))
                ; cropping
                binary-signal (net.imglib2.view.Views/offset (net.imglib2.view.Views/interval binary-signal offset stop-point)
                                                             (long-array start-point))
                cell-binary-signals (net.imglib2.view.Views/offset (net.imglib2.view.Views/interval cell-binary-signals offset stop-point)
                                                                   (long-array start-point))

                ;inverted (normalize-in-mask inverted mask)
                crop-name (str "crop_sX_" start-x "_crop_sY_" start-y)

                _ (println (str basename "_" crop-name "_filtered.tif"))
                filtered-regions (ij/open-img (str basename "_" crop-name "_filtered.tif"))

                crop-imgs (xycz-to-xyz filtered-regions)]

            ; Calculating stats
            ; Signal stats
            (img/map-img! (fn [label signal]
                            (let [lid (str start-x "_" start-y "_" (.get (.get label)))]
                              (reset! binary-signal-counts
                                      (assoc @binary-signal-counts
                                        lid (if (get @binary-signal-counts lid)
                                              (+ (get @binary-signal-counts lid) (.get (.get signal)))
                                              (.get (.get signal)))))))
                          (first crop-imgs) binary-signal))))
      ;; binary signal render
      (doall
        (for [start-x (range 0 (first img-dims) (first crop-dims))
              start-y (range 0 (second img-dims) (second crop-dims))]
          (let [start-point [start-x start-y 0]
                _ (println "Starting crop " start-point)
                offset (long-array start-point)
                stop-point (long-array (map #(dec (min (+ %1 %2)
                                                       %3))
                                            offset
                                            crop-dims
                                            img-dims))

                cell-binary-signals (net.imglib2.view.Views/offset (net.imglib2.view.Views/interval cell-binary-signals offset stop-point)
                                                                   (long-array start-point))

                ;inverted (normalize-in-mask inverted mask)
                crop-name (str "crop_sX_" start-x "_crop_sY_" start-y)

                _ (println (str basename "_" crop-name "_filtered.tif"))
                filtered-regions (ij/open-img (str basename "_" crop-name "_filtered.tif"))

                crop-imgs (xycz-to-xyz filtered-regions)]

            ; Calculating stats
            ; Signal stats
            (img/map-img! (fn [label signal]
                            (let [lid (str start-x "_" start-y "_" (.get (.get label)))]
                              (.set (.get signal) (double (get @binary-signal-counts lid)))))
                          (first crop-imgs) cell-binary-signals)))))
    #_(let [ks (keys @voxel-counts)]
      (spit (str basename "_cellStats.csv")
            (with-out-str
              (println (string/join "," ["cellID" "x" "y" "z" "voxels" "signalSum" "min_x" "min_y" "min_z" "max_x" "max_y" "max_z" "binarySignal"]))
              (doseq [k ks]
                (when (> (get @crop-id k) 0.001)
                  (println (string/join "," [k
                                             (float (/ (get @cell-xs k) (get @voxel-counts k)))
                                             (float (/ (get @cell-ys k) (get @voxel-counts k)))
                                             (float (/ (get @cell-zs k) (get @voxel-counts k)))
                                             (get @voxel-counts k) (get @signal-counts k)
                                             (get @min-xs k) (get @min-ys k) (get @min-zs k)
                                             (get @max-xs k) (get @max-ys k) (get @max-zs k)
                                             (get @binary-signal-counts k)])))))))
    (ij/save-img (fun.imagej.ops.convert/float32 membrane) (str basename "_membrane.tif"))
    (ij/save-img (fun.imagej.ops.convert/float32 signal) (str basename "_signal.tif"))
    (ij/save-img (fun.imagej.ops.convert/float32 cell-labels) (str basename "_cellLabels.tif"))
    (ij/save-img (fun.imagej.ops.convert/float32 cell-signals) (str basename "_cellSignals.tif"))
    (ij/save-img (fun.imagej.ops.convert/float32 cell-binary-signals) (str basename "_cellBinarySignals.tif"))
    (ij/save-img (fun.imagej.ops.convert/float32 cell-counts) (str basename "_cellCounts.tif"))))

;(-main)

;(def inverted (-main))
;(def rcomp (region-competition inverted))

;(in-ns 'incisor-cell-counter.full-analysis)

#_(def mask
  (let [full-mask (fun.imagej.ops.threshold/intermodes
                    (img/concat-imgs
                      (doall (for [k (range (img/get-size-dimension c1 2))]
                               (let [source (img/hyperslice c1 2 k)
                                     output (fun.imagej.ops.create/img source)]
                                 (fun.imagej.ops.filter/gauss output source (double-array [5 5]))
                                 output
                                 #_(fun.imagej.ops.threshold/triangle output))))))]
    (fun.imagej.ops.morphology/erode full-mask (shape/sphere 3))))
#_(ij/show mask "mask")

#_(def masked-signal
  (ij/show
    (fun.imagej.ops.math/multiply (fun.imagej.ops.convert/uint16 mask)
                                  (fun.imagej.ops.convert/uint16 (fun.imagej.ops.threshold/maxEntropy c2)))))
