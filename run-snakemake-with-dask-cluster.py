import subprocess
import sys

import dask
import yaml

if __name__ == "__main__":
    with open("config.yaml") as f:
        config = yaml.safe_load(f)

    for key, value in config["dask"].items():
        dask.config.set({key: value})

    # Make a dask cluster from the module.class string specified in config
    cluster_class = config["cluster"]["class"]
    cluster_module, cluster_class = cluster_class.rsplit(".", 1)
    cluster_module = __import__(cluster_module, fromlist=[cluster_class])
    cluster_class = getattr(cluster_module, cluster_class)
    cluster = cluster_class(**config["cluster"].get("kwargs", {}))
    print(
        f"Created cluster {cluster.__class__.__name__} at {cluster.scheduler_address}"
    )
    print(f"Dashboard at {cluster.dashboard_link}")
    if "adapt" in config["cluster"]:
        cluster.adapt(**config["cluster"]["adapt"])
        print(f"Scaling cluster with {config['cluster']['adapt']}")

    # Run snakemake with the args that this script was called with, adding in the
    # scheduler address as a config we need to run it in a subprocess, then when that
    # is done shutdown the dask cluster

    subprocess.run(
        [
            "snakemake",
        ]
        + sys.argv[1:]
        + ["--config", f"scheduler_address={cluster.scheduler_address}"],
        check=True,
    )
    print(f"Shutting down cluster. {len(cluster.workers)} workers live.")
    cluster.close()
