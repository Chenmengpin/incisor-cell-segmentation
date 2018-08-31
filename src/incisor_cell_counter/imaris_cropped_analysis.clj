(ns incisor-cell-counter.imaris-cropped-analysis
  (:gen-class)
  (:import (mosaic.region_competition.RC AlgorithmRC)
           (mosaic.region_competition.DRS SobelVolume)
           (ij.process StackStatistics)
           (mosaic.core.imageUtils.images IntensityImage)
           (mosaic.plugins Region_Competition)
           (net.imglib2.view Views)
           (loci.plugins.in ImporterOptions)
           (net.imglib2.img.array ArrayImgs))
  (:require [fun.imagej.img :as img]
            [fun.imagej.imp :as ij1]
            [fun.imagej.conversion :as convert]
            [fun.imagej.img.cursor :as cursor]
            [fun.imagej.img.shape :as shape]
            [fun.imagej.img.type :as tpe]
            [fun.imagej.core :as ij]
            [flatland.ordered.set :as os]
    ;[fun.imagej.ops :as ops]
            [fun.imagej.img.utils :as img-utils]
            [fun.imagej.mesh :as msh]
            [fun.imagej.segmentation.fast-segmentation :as fseg]
            [clojure.java.io :as io]
            [clojure.string :as string]))

;(def display false)
(def display false)

#_(when-not display
  (.setProperty (System/getProperties)
                "java.awt.headless" "true"))

(when display
  (defonce ijui (ij/show-ui)))

#_(defn normal-form-labels
  "Relabel starting at 1, 0 is treated separately, pass through, as background"
  [label-image]
  (let [ordered-labels (atom (os/ordered-set))
        label-img (convert/imp->img (ArrayImgs/ints (.getDataLabel label-image) (long-array (seq (img/dimensions input)))))]
    (img/map-img! (fn [cur]
                    (swap! ordered-labels conj (.get (.get cur))))
                  label-img)
    (let [label-map (zipmap @ordered-labels (range))]
      #_(first (img/map-img! (fn [dest source]
                               (.set (.get dest)
                                     (get label-map (.get (.get source)))))
                             (fun.imagej.ops.create/img label-img) label-img))
      (first
        label-img))))

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

(defn region-competition
    "Run region competition on an input image."
    [input output-prefix voxel-res]
    (let [input-imp (convert/img->imp (fun.imagej.ops.convert/float32 input))
          cal (ij1/get-calibration input-imp)
          ; Imaris RL2 (hard coded to Amnon's first zstep)
          ;_ (set! (.pixelWidth cal) 0.9848675)
          ;_ (set! (.pixelHeight cal) 0.9848675)
          ;_ (set! (.pixelDepth cal) 1.9991803)
          _ (set! (.pixelWidth cal) (:w voxel-res))
          _ (set! (.pixelHeight cal) (:h voxel-res))
          _ (set! (.pixelDepth cal) (:d voxel-res))
          intensity-image (mosaic.core.imageUtils.images.IntensityImage. input-imp false)
          label-image (mosaic.core.imageUtils.images.LabelImage. (fun.imagej.ops.convert/int16 (fun.imagej.ops.create/img input)))
          ;label-initializer (mosaic.region_competition.initializers.BubbleInitializer. label-image); use MaximaBubbles
          ;_ (.initialize label-initializer 4 6)
          label-initializer (mosaic.region_competition.initializers.MaximaBubbles. intensity-image label-image
                                                                                   0.5 ;2.0
                                                                                   0.00005 ;0.005
                                                                                   2
                                                                                   2); use MaximaBubbles
          _ (.initialize label-initializer)
          settings (mosaic.region_competition.RC.Settings.)
          external-energy (mosaic.region_competition.energies.E_PS. label-image intensity-image 5 0.0001)
          internal-energy (mosaic.region_competition.energies.E_CurvatureFlow. label-image 3 cal)
          merge-energy (mosaic.region_competition.energies.E_KLMergingCriterion. 0 0.005)
          image-model (mosaic.region_competition.energies.ImageModel. external-energy
                                                                      internal-energy
                                                                      merge-energy
                                                                      settings)
          ;algorithm-rc (mosaic.region_competition.RC.AlgorithmRC. intensity-image label-image image-model settings)
          label-array (.getDataLabel label-image)
          display-label (ArrayImgs/ints label-array (long-array (seq (img/dimensions input))))]
      (if (zero? (reduce + label-array))
        display-label
        (do
                ;display-label (.convertToImg label-image "label progress")]
            (when display (ij/show input "input"))
            (when display
              (ij/show display-label))

                #_(ij1/save-imp-as-tiff display-label
                                        (str output-prefix "_" 0 ".tif"))
            #_(ij/save-img display-label
                           (str output-prefix "_" 0 ".tif"))
                #_(ij/save-img (fun.imagej.ops.convert/uint16 input)
                             (str output-prefix "_input.tif"))
                ;; Previous working
                #_(dotimes [num-step 16];100
                    (println "Step" num-step)
                    ;(.performIteration algorithm-rc)
                    (.performIteration (mosaic.region_competition.RC.AlgorithmRC. intensity-image label-image image-model settings))
                    ;(normal-form-labels label-image)
                    ; Relabel label-image
                    ;(when display (ij1/update-imp! display-label))
                    #_(ij1/save-imp-as-tiff display-label
                                            (str output-prefix "_" num-step ".tif"))
                    #_(ij/save-img display-label
                                   (str output-prefix "_" num-step ".tif")))
                (loop [num-step 50
                       convergence? false];100
                  (println "Step" num-step)
                  (when-not (or convergence?
                                (neg? num-step))
                    (recur (dec num-step)
                           (.performIteration (mosaic.region_competition.RC.AlgorithmRC. intensity-image label-image image-model settings)))))
            display-label))))

; For debugging
;(def display false)
#_(let [filtered-regions #_(ij/open-img "/Users/kharrington/Data/Sharir_Amnon/Results008/24h-244c_crop_sX_300_crop_sY_200_filtered.tif")
        (ij/open-img "/Users/kharrington/Data/Sharir_Amnon/Results011/45m-0062c_crop_sX_100_crop_sY_100_filtered.tif")
        crop-imgs (xycz-to-xyz filtered-regions)
        membrane (nth crop-imgs 3)]
    (region-competition membrane "/Users/kharrington/Data/Sharir_Amnon/Results011/45m-0062c_crop_sX_100_crop_sY_100_test" {:w 1 :h 1 :d 2}))

; no maxima: 45m-007c_crop_sX_500_crop_sY_500_filtered.tif

; For debugging from final crops
#_(do
    (ij/show-ui)
    (def filtered-regions (ij/open-img "/Users/kharrington/Data/Sharir_Amnon/Results013/45m-0062d_crop_sX_100_crop_sY_100_filtered.tif"))
    (ij/show filtered-regions)
    (def membrane (img/hyperslice filtered-regions 2 3))
    (def signal (img/hyperslice filtered-regions 2 4))
    (def output-prefix "/Users/kharrington/Data/Sharir_Amnon/Results013/45m-0062d_crop_sX_100_crop_sY_100_debug")
    (def voxel-res {:w 1 :h 1 :d 2})
    (def mask (let [full-mask (fun.imagej.ops.threshold/intermodes
                                (img/concat-imgs
                                  (doall (for [k (range (img/get-size-dimension membrane 2))]
                                           (let [source (img/hyperslice membrane 2 k)
                                                 output (fun.imagej.ops.create/img source)]
                                             (fun.imagej.ops.filter/gauss output source (double-array [5 5]))
                                             output)))))]
                (fun.imagej.ops.morphology/erode full-mask (shape/sphere 5))))
    (def inverted (fun.imagej.ops.math/multiply (fun.imagej.ops.convert/uint16 (fun.imagej.ops.image/invert (fun.imagej.ops.create/img membrane) membrane))
                                                (fun.imagej.ops.convert/uint16 mask)))
    (def regions (region-competition inverted output-prefix voxel-res))
    )

; For debugging at whole image scale; testing externally provided maxima
#_(do
  (ij/show-ui)
  (def membrane (ij/open-img "/Users/kharrington/Data/Sharir_Amnon/Results013/45m-0062d_membrane.tif"))
  (ij/show membrane)
  (def signal (ij/open-img "/Users/kharrington/Data/Sharir_Amnon/Results013/45m-0062d_signal.tif"))
  (def output-prefix "/Users/kharrington/Data/Sharir_Amnon/Results013/45m-0062d_debug")
  (def voxel-res {:w 1 :h 1 :d 2})
  (def mask (let [full-mask (fun.imagej.ops.threshold/triangle
                              (img/concat-imgs
                                (doall (for [k (range (img/get-size-dimension membrane 2))]
                                         (let [source (img/hyperslice membrane 2 k)
                                               output (fun.imagej.ops.create/img source)]
                                           (fun.imagej.ops.filter/gauss output source (double-array [5 5]))
                                           output)))))]
              (fun.imagej.ops.morphology/erode full-mask (shape/sphere 5))))
  (ij/show mask)
  (def inverted (fun.imagej.ops.math/multiply (fun.imagej.ops.convert/uint16 (fun.imagej.ops.image/invert (fun.imagej.ops.create/img membrane) membrane))
                                              (fun.imagej.ops.convert/uint16 mask)))
  (def start-x 100)
  (def start-y 100)
  (def img-dims (seq (img/dimensions signal)))
  (def crop-dims [100 100 (last img-dims)])
  (def region (process-cropped-region start-x start-y crop-dims img-dims output-prefix  membrane signal mask inverted voxel-res))
  (ij/save-img region (str output-prefix "_regions_002.tif"))
  (def cc-region (color-code-by-size region))
  (ij/save-img cc-region (str output-prefix "_regions_003.tif"))
  )


(defn filter-regions
  "Filter regions"
  [label-image counts-img signal-img]
  (let [;label-img (convert/imp->img label-image)
        label-img label-image
        output-labels (img/create-img-like label-img)
        output-signal (img/create-img-like label-img) #_(fun.imagej.ops.convert/int16 (fun.imagej.ops.create/img signal-img))
        signal-table (atom {})
        min-size Double/MIN_VALUE #_70
        max-size Double/MAX_VALUE #_300]
    (println "out-signal" (.firstElement output-signal) (.firstElement label-img))
    (img/map-img! (fn [in-label in-signal]
                    (swap! signal-table update (.get (.get in-label)) #(if %
                                                                         (+ % (.get (.get in-signal)))
                                                                         (.get (.get in-signal)))))
                  label-img signal-img); Compute the signal table
    (img/map-img! (fn [in-label in-signal in-counts out-label out-signal]
                    (let [c (.get (.get in-counts))]
                      (when (and (< c max-size)
                               (> c min-size))
                        (.set (.get out-label) (.get (.get in-label)))
                        (.set (.get out-signal) (int (get @signal-table
                                                      (.get (.get in-label))))))))
                  label-img signal-img counts-img output-labels output-signal)
    [output-labels output-signal counts-img]))

#_(defn region-competition
    "Run region competition on an input image."
    [input]
    (let [input-imp (convert/img->imp input)
          cal (ij1/get-calibration input-imp)
          intensity-image (mosaic.core.imageUtils.images.IntensityImage. input-imp true)
          label-image (mosaic.core.imageUtils.images.LabelImage. (fun.imagej.ops.create/img input))
          label-initializer (mosaic.region_competition.initializers.MaximaBubbles. intensity-image label-image
                                                                                   0.25 ;2.0
                                                                                   0.0005 ;0.005
                                                                                   3
                                                                                   4); use MaximaBubbles
          _ (.initialize label-initializer)
          settings (mosaic.region_competition.RC.Settings.)
          external-energy (mosaic.region_competition.energies.E_PS. label-image intensity-image 14 0.04)
          internal-energy (mosaic.region_competition.energies.E_CurvatureFlow. label-image 8 cal)
          merge-energy (mosaic.region_competition.energies.E_KLMergingCriterion. 0 0.02)
          image-model (mosaic.region_competition.energies.ImageModel. external-energy
                                                                      internal-energy
                                                                      merge-energy
                                                                      settings)
          sobel-imp (ij.ImagePlus. "sobelInput" (.convertToFloat (.getImageStack (.convertToImg intensity-image "intensity"))))
          sobel-vol (SobelVolume. sobel-imp)
          _ (.sobel3D sobel-vol)
          sobel-ip (ij.ImagePlus. "xxxx" (.getImageStack sobel-vol))
          ;ss (StackStatistics. sobel-ip)
          ;so displayrange
          edge-image (IntensityImage. sobel-ip)
          _ (println intensity-image label-image edge-image)
          algorithm-rc (mosaic.region_competition.DRS.AlgorithmDRS. intensity-image label-image edge-image image-model (mosaic.region_competition.DRS.Settings.))]
      ;algorithm-rc (mosaic.region_competition.RC.AlgorithmRC. intensity-image label-image image-model settings)]
      (ij/show input "input")
      (ij1/show-imp (.convertToImg label-image "label progress"))
      (dotimes [num-step 10]
        (println "Step" num-step)
        (.performIteration algorithm-rc))
      label-image))

(defn color-code-by-size
  "Color code a label image by size"
  [input]
  (let [region-counts (frequencies (img/img-to-seq input))]
    (second (img/map-img (fn [icur ocur]
                          (.set (.get ocur)
                                (get region-counts
                                     (.get (.get icur)))))
                        input (fun.imagej.ops.create/img input)))
    ; Compute size lookup-table
    ; Make image that stores count of lookup table
    ))

(defn normalize-in-mask
  "BROKEN/INCOMPLETE"
  [input mask]
  (let [minv (atom Double/MAX_VALUE)
        maxv (atom Double/MIN_VALUE)
        totmin (atom Double/MAX_VALUE)]
    (first (img/map-img (fn [icur mcur]
                          (if (.get (.get mcur))
                            (do (reset! minv (min @minv (.get (.get icur))))
                                (reset! maxv (min @maxv (.get (.get icur)))))
                            (reset! totmin (min @totmin (.get (.get icur))))))
                        input mask))))

#_(defn interleave-imgs
  "Interleave a seq o images"
  [imgs]
    (let [dimensions (dimensions (first imgs))
          stack (fun.imagej.ops.create/img (long-array (concat dimensions
                                                               [(count imgs)])))]
      (dotimes [k (count imgs)]
        (map-img cursor/copy-real
                 (Views/hyperSlice stack (count dimensions) k)
                 (nth imgs k)))
      stack))

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

(defn process-cropped-region
  [start-x start-y crop-dims img-dims basename membrane signal mask inverted voxel-res]
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
        _ (println (seq (img/dimensions membrane)) start-point (seq stop-point) img-dims)
        signal (net.imglib2.view.Views/offset (net.imglib2.view.Views/interval signal offset stop-point)
                                              (long-array start-point))
        mask (net.imglib2.view.Views/offset (net.imglib2.view.Views/interval mask offset stop-point)
                                            (long-array start-point))
        inverted (net.imglib2.view.Views/offset (net.imglib2.view.Views/interval inverted offset stop-point)
                                                (long-array start-point))
        ;inverted (normalize-in-mask inverted mask)
        crop-name (str "crop_sX_" start-x "_crop_sY_" start-y)
        rcomp (try (region-competition inverted (str basename "_" crop-name "_label") voxel-res)
                   (catch Exception e (do (println "ERROR: " basename "_" crop-name " region competition failure")
                                          (.printStackTrace e)
                                          (fun.imagej.ops.create/img inverted))))
        counts (color-code-by-size rcomp)
        filtered-regions (fun.imagej.ops.convert/float32
                           (xyz-to-xycz (conj (filter-regions rcomp counts signal)
                                              membrane signal)))]
    ;(ij/save-img filtered-regions (str basename "_" crop-name "_filtered.tif"))
    (ij1/save-imp-as-tiff (convert/img->imp filtered-regions)
                          (str basename "_" crop-name "_filtered.tif"))
    #_(when display (ij1/show-imp filtered-regions))
    (when display (ij/show filtered-regions))
    #_(when display (ij1/show-imp rcomp))
    (when display (ij/show rcomp))
    (ij/save-img (first (img/map-img cursor/copy
                                     (fun.imagej.ops.create/img inverted)
                                     inverted))
                 (str basename "_" crop-name "_inverted.tif"))
    (ij/save-img counts
                 (str basename "_" crop-name "_counts.tif"))
    (ij/save-img rcomp
                 (str basename "_" crop-name "_label.tif"))))

(defn -main
  [& args]
  (println "Starting")
  (let [directory (or (first args)
                      "/Users/kharrington/Data/Sharir_Amnon/Results013/")
        just-name (or (second args)
                      "45m-0062d")
        basename (str directory just-name)
        ;c1-filename (str directory "C1-" just-name ".tif")
        ;c2-filename (str directory "C2-" just-name ".tif")
        ;membrane (ij/open-img c1-filename)
        ;signal (ij/open-img c2-filename)
        imps (seq (load-imaris-rl2 (str directory just-name ".ims")))
        voxel-res (let [cal (ij1/get-calibration (first imps))]
                    {:w (.pixelWidth cal)
                     :h (.pixelHeight cal)
                     :d (.pixelDepth cal)})
        _ (println voxel-res)
        [membrane signal] (map convert/imp->img imps)
        ;_ (ij1/show-imp membrane)
        ;_ (ij1/show-imp signal)
        _ (println "Images opened")
        ;mask (let [full-mask (fun.imagej.ops.threshold/intermodes
        ;                       (img/concat-imgs
        ;                         (doall (for [k (range (img/get-size-dimension membrane 2))]
        ;                                  (let [source (img/hyperslice membrane 2 k)
        ;                                        output (fun.imagej.ops.create/img source)]
        ;                                    (fun.imagej.ops.filter/gauss output source (double-array [5 5]))
        ;                                    output
        ;                                    #_(fun.imagej.ops.threshold/triangle output))))))]
        ;       (fun.imagej.ops.morphology/erode full-mask (shape/sphere 3)))

        ;; Cropping
        img-dims (seq (img/dimensions membrane))
        ;crop-dims [100 100 50]
        ;crop-dims [200 200 50]
        #_start-point #_(map #(+ %1 (rand-int (- %2 (* 2 %1))))
                         crop-dims img-dims)
        mask (let [full-mask (fun.imagej.ops.threshold/triangle ;previously intermodes
                               (img/concat-imgs
                                 (doall (for [k (range (img/get-size-dimension membrane 2))]
                                          (let [source (img/hyperslice membrane 2 k)
                                                output (fun.imagej.ops.create/img source)]
                                            (fun.imagej.ops.filter/gauss output source (double-array [5 5]))
                                            output
                                            #_(fun.imagej.ops.threshold/triangle output))))))]
               (fun.imagej.ops.morphology/erode
                 full-mask (shape/sphere 5)))
        _ (println "Mask")
        inverted (fun.imagej.ops.math/multiply (fun.imagej.ops.image/invert (fun.imagej.ops.create/img membrane) membrane)
                                               (fun.imagej.ops.convert/uint16 mask))
        histogram (fun.imagej.ops.image/histogram inverted)
        min-val (.createVariable (.firstElement inverted))
        _ (.getCenterValue histogram 1 min-val)
        inverted (fun.imagej.ops.math/multiply (fun.imagej.ops.convert/uint16 mask)
                                               (fun.imagej.ops.math/subtract inverted min-val))
        _ (println "Starting region competition")
        ;crop-dims [30 30 (last img-dims)]]
        ;crop-dims [50 50 (last img-dims)]
        crop-dims [100 100 (last img-dims)]
        ]
    (doall
      ; For processing extra crops
      ;(for [[start-x start-y start-z] (map #(map read-string (string/split % #" "))
      ;                                     (string/split-lines (slurp "/mnt/lfs2/kharrington/incisor/imaris_cropped/242_crop_errors.clj")))]
      ; For proper runs
      (for [start-x (range 0 (first img-dims) (first crop-dims))
            start-y (range 0 (second img-dims) (second crop-dims))]
      ; For testing
      ;(for [start-x (range 100 200 (first crop-dims))
      ;      start-y (range 100 200 (second crop-dims))]
      ; Encapsulated version of prev working
      (process-cropped-region start-x start-y crop-dims img-dims basename  membrane signal mask inverted voxel-res)
      ; Previously working
      #_(let [start-point [start-x start-y 0]
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
            _ (println (seq (img/dimensions membrane)) start-point (seq stop-point) img-dims)
            signal (net.imglib2.view.Views/offset (net.imglib2.view.Views/interval signal offset stop-point)
                                                  (long-array start-point))
            mask (net.imglib2.view.Views/offset (net.imglib2.view.Views/interval mask offset stop-point)
                                                (long-array start-point))
            inverted (net.imglib2.view.Views/offset (net.imglib2.view.Views/interval inverted offset stop-point)
                                                    (long-array start-point))
            ;inverted (normalize-in-mask inverted mask)
            crop-name (str "crop_sX_" start-x "_crop_sY_" start-y)
            rcomp (try (region-competition inverted (str basename "_" crop-name "_label") voxel-res)
                       (catch Exception e (do (println "ERROR: " basename "_" crop-name " region competition failure")
                                              (.printStackTrace e)
                                              (fun.imagej.ops.create/img inverted))))
            counts (color-code-by-size rcomp)
            filtered-regions (fun.imagej.ops.convert/float32
                               (xyz-to-xycz (conj (filter-regions rcomp counts signal)
                                                  membrane signal)))]
        ;(ij/save-img filtered-regions (str basename "_" crop-name "_filtered.tif"))
        (ij1/save-imp-as-tiff (convert/img->imp filtered-regions)
                              (str basename "_" crop-name "_filtered.tif"))
        #_(when display (ij1/show-imp filtered-regions))
        (when display (ij/show filtered-regions))
        #_(when display (ij1/show-imp rcomp))
        (when display (ij/show rcomp))
        (ij/save-img (first (img/map-img cursor/copy
                                         (fun.imagej.ops.create/img inverted)
                                         inverted))
                     (str basename "_" crop-name "_inverted.tif"))
        (ij/save-img counts
                     (str basename "_" crop-name "_counts.tif"))
        (ij/save-img rcomp
                     (str basename "_" crop-name "_label.tif")))))))

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
