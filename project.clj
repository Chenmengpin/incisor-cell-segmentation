(defproject incisor-cell-counter "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :jvm-opts ["-Xmx25g" "-server"]
  :main incisor-cell-counter.imaris-cropped-analysis
  :dependencies [[org.clojure/clojure "1.8.0"]
                 [fun.imagej/fun.imagej "0.2.5"]
                 [mosaic/MosaicSuite "1.0.8_Full"]
                 [org.flatland/ordered "1.5.6"]
                 ]
  :repositories [["imagej" "http://maven.imagej.net/content/groups/hosted/"]
                 ["imagej-releases" "http://maven.imagej.net/content/repositories/releases/"]
                 ["ome maven" "http://artifacts.openmicroscopy.org/artifactory/maven/"]
                 ["imagej-snapshots" "http://maven.imagej.net/content/repositories/snapshots/"]
                 ["sonatype-snapshots" "https://oss.sonatype.org/content/repositories/snapshots/"]
                 ["snapshots" {:url "https://clojars.org/repo"
                               :username :env/CI_DEPLOY_USERNAME
                               :password :env/CI_DEPLOY_PASSWORD
                               :sign-releases false}]
                 ["releases" {:url "https://clojars.org/repo"
                              :username :env/CI_DEPLOY_USERNAME
                              :password :env/CI_DEPLOY_PASSWORD
                              :sign-releases false}]
                 ])
(require 'cemerick.pomegranate.aether)
(cemerick.pomegranate.aether/register-wagon-factory!
  "http" #(org.apache.maven.wagon.providers.http.HttpWagon.))
