import os
import shutil
import pandas as pd

path = "results/daily/{}_{}.json"
template = "snakemake -call{} --configfile config/config.yaml -- {}"

target = "results/live/now.json"

if __name__ == "__main__":

    now = pd.Timestamp.now()
    day = now.strftime("%Y-%m-%d")
    period = now.hour * 2 + now.minute // 30

    outfile = path.format(day, period)

    os.system(template.format(" --touch", outfile))
    os.system(template.format("", outfile))

    shutil.copy(outfile, "results/live/latest.json")
    os.remove(outfile)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(10, 3.5))

    ax.text(0.5, 0.5, f"{day} {period}", fontsize=12, ha='center')

    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xticks([])
    ax.set_yticks([])

    plt.savefig('results/live/current_period.png')
    plt.show()
