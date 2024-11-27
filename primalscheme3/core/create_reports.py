import pathlib

import plotly.graph_objects as go
from plotly.subplots import make_subplots


def generate_all_plots_html(
    plot_data: dict, outdir: pathlib.Path, offline_plots: bool = True
) -> str:
    """Generate all the plots for a scheme from the plot_data"""
    # Generate the plot for each MSA
    plot_html = []

    for index, (chromname, data) in enumerate(plot_data.items()):
        plot_html.append(
            generate_plot_html(
                chromname=chromname,
                msa_data=data,
                outdir=outdir,
                offline_plots=True if offline_plots and index == 0 else False,
            )
        )
    return "\n".join(plot_html)


def generate_plot_html(
    chromname: str, msa_data: dict, outdir: pathlib.Path, offline_plots=True
) -> str:
    # Create an empty figure with the Fprimer hovers
    fig = make_subplots(
        cols=1,
        rows=4,
        shared_xaxes=True,
        vertical_spacing=0.02,
        row_heights=[0.2, 0.2, 0.2, 0.1],
        specs=[
            [{"secondary_y": False}],
            [{"secondary_y": True}],
            [{"secondary_y": False}],
            [{"secondary_y": False}],
        ],
    )

    # Extract amplicon data from the msa_data
    # Filter primers that are circular
    circular_pp = []
    amplicons = []
    for _, pp in msa_data["amplicons"].items():
        if pp["cs"] > pp["ce"]:
            circular_pp.append(pp)
        else:
            amplicons.append(pp)

    length = msa_data["dims"][1]

    if len(amplicons) == 0:
        npools = 0
    else:
        npools = max([x["p"] for x in amplicons])

    fig.add_trace(
        go.Scattergl(
            x=[x["cs"] for x in amplicons],
            y=[x["p"] for x in amplicons],
            opacity=0,
            name="FPrimers",
            hovertext=[f"{x['n']}_LEFT" for x in amplicons],
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scattergl(
            x=[x["ce"] for x in amplicons],
            y=[x["p"] for x in amplicons],
            opacity=0,
            name="RPrimers",
            hovertext=[f"{x['n']}_RIGHT" for x in amplicons],
        ),
        row=1,
        col=1,
    )
    # Add the uncovered regions
    for start, stop in msa_data["uncovered"].items():
        fig.add_vrect(
            x0=start,
            x1=stop + 1,
            fillcolor="#F0605D",
            line=dict(width=0),
            opacity=0.5,
            row=1,  # type: ignore
            col=1,  # type: ignore
        )

    # Plot the amplicons lines
    for amplicon in amplicons:
        fig.add_shape(
            type="line",
            y0=amplicon["p"],
            y1=amplicon["p"],
            x0=amplicon["cs"],
            x1=amplicon["ce"],
            line=dict(color="LightSeaGreen", width=5),
            row=1,
            col=1,
        )
        fig.add_shape(
            type="rect",
            y0=amplicon["p"] - 0.05,
            y1=amplicon["p"] + 0.05,
            x0=amplicon["s"],
            x1=amplicon["cs"],
            fillcolor="LightSalmon",
            line=dict(color="darksalmon", width=2),
            row=1,
            col=1,
        )
        fig.add_shape(
            type="rect",
            y0=amplicon["p"] - 0.05,
            y1=amplicon["p"] + 0.05,
            x0=amplicon["ce"],
            x1=amplicon["e"],
            fillcolor="LightSalmon",
            line=dict(color="darksalmon", width=2),
            row=1,
            col=1,
        )

    # Plot the circular primers
    for pp in circular_pp:
        # Add the left side line
        fig.add_shape(
            type="line",
            y0=pp["p"],
            y1=pp["p"],
            x1=pp["ce"],
            x0=0,
            line=dict(color="LightSeaGreen", width=5),
            row=1,
            col=1,
        )
        # Add the right side line
        fig.add_shape(
            type="line",
            y0=pp["p"],
            y1=pp["p"],
            x0=pp["cs"],
            x1=length,
            line=dict(color="LightSeaGreen", width=5),
            row=1,
            col=1,
        )
        fig.add_shape(
            type="rect",
            y0=pp["p"] - 0.05,
            y1=pp["p"] + 0.05,
            x0=pp["s"],
            x1=pp["cs"],
            fillcolor="LightSalmon",
            line=dict(color="LightSalmon", width=2),
            row=1,
            col=1,
        )
        fig.add_shape(
            type="rect",
            y0=pp["p"] - 0.05,
            y1=pp["p"] + 0.05,
            x0=pp["ce"],
            x1=pp["e"],
            fillcolor="LightSalmon",
            line=dict(color="LightSalmon", width=2),
            row=1,
            col=1,
        )

    # Add the base occupancy
    occupancy_data = [
        (int(index), float(oc)) for index, oc in msa_data["occupancy"].items()
    ]
    fig.add_trace(
        go.Scattergl(
            x=[x[0] for x in occupancy_data],
            y=[x[1] for x in occupancy_data],
            mode="lines",
            name="Base Occupancy",
            line=dict(color="#F0605D", width=2),
            fill="tozeroy",
            opacity=0.5,
        ),
        row=2,
        col=1,
    )

    ## Plot the GC data
    gc_data = [(int(index), float(gc)) for index, gc in msa_data["gc"].items()]
    fig.add_trace(
        go.Scattergl(
            x=[x[0] for x in gc_data],
            y=[x[1] for x in gc_data],
            mode="lines",
            name="GC Prop",
            line=dict(color="#005c68", width=2),
        ),
        row=2,
        col=1,
        secondary_y=True,
    )

    # Add the entropy plot
    entropy_data = [
        (int(index), float(entropy)) for index, entropy in msa_data["entropy"].items()
    ]
    fig.add_trace(
        go.Scattergl(
            x=[x[0] for x in entropy_data],
            y=[x[1] for x in entropy_data],
            opacity=1,
            name="Sequence Entropy",
            mode="lines",
        ),
        row=3,
        col=1,
    )

    # Add all possible Fkmers
    fkmer_data = [
        (end, num_seqs) for end, num_seqs in msa_data["thermo_pass"]["F"].items()
    ]
    fig.add_trace(
        go.Scattergl(
            x=[x[0] for x in fkmer_data],
            y=[1 for _ in fkmer_data],
            hovertext=[x[1] for x in fkmer_data],
            marker=dict(symbol="triangle-right", size=10),
            mode="markers",
            name="Passing Forward Primers",
            hovertemplate="Number Seqs: %{hovertext}",
        ),
        row=4,
        col=1,
    )
    # Add all possible Rkmers
    rkmer_data = [
        (start, num_seqs) for start, num_seqs in msa_data["thermo_pass"]["R"].items()
    ]
    fig.add_trace(
        go.Scattergl(
            x=[x[0] for x in rkmer_data],
            y=[0.5 for _ in rkmer_data],
            hovertext=[x[1] for x in rkmer_data],
            marker=dict(symbol="triangle-left", size=10),
            mode="markers",
            name="Passing Reverse Primers",
            hovertemplate="Number Seqs: %{hovertext}",
        ),
        row=4,
        col=1,
    )

    # Add the regions
    # If regions are present add them
    if "regions" in msa_data:
        data = msa_data.get("regions", [])
        for region in data:
            fig.add_shape(
                type="line",
                y0=1.5,
                y1=1.5,
                x0=region["s"],
                x1=region["e"],
                fillcolor="Green",
                line=dict(color="Green", width=5),
                row=1,
                col=1,
            )
        fig.add_trace(
            go.Scatter(  # doesn't need WebGL
                x=[x["s"] + ((x["e"] - x["s"]) // 2) for x in data],
                y=[1.5 for _ in data],
                name="Regions",
                mode="markers",
                text=[x["n"] for x in data],
                hovertemplate="%{text}",
                opacity=0,  # Make the markers invisible
            ),
        )

    # Add the base plot settings
    fig.update_xaxes(
        showline=True,
        mirror=True,
        ticks="outside",
        linewidth=2,
        linecolor="black",
        tickformat=",d",
        title_font=dict(size=18, family="Arial", color="Black"),
        range=[0, length],
        title="",  # Blank title for all x-axes
    )
    fig.update_yaxes(
        showline=True,
        mirror=True,
        ticks="outside",
        linewidth=2,
        linecolor="black",
        fixedrange=True,
        title_font=dict(size=18, family="Arial", color="Black"),
    )

    # Update the top plot
    fig.update_yaxes(
        range=[0.5, npools + 0.5],
        title="pool",
        tickmode="array",
        tickvals=sorted({x["p"] for x in amplicons}),
        row=1,
        col=1,
    )
    # Update the second plot
    fig.update_yaxes(
        range=[-0.1, 1.1],
        title="Base Occupancy",
        tickmode="array",
        tickvals=[0, 0.25, 0.5, 0.75, 1],
        row=2,
        col=1,
        secondary_y=False,
    )
    fig.update_yaxes(
        range=[-0.1, 1.1],
        title="GC%",
        tickmode="array",
        tickvals=[0, 0.25, 0.5, 0.75, 1],
        ticktext=[0, 25, 50, 75, 100],
        row=2,
        col=1,
        secondary_y=True,
        side="right",
    )
    # Update the third plot
    fig.update_yaxes(
        title="Entropy",
        tickmode="array",
        row=3,
        col=1,
    )
    # Update the fourth plot
    fig.update_yaxes(
        title="Thermo-passing Primers",
        range=[0.5 - 0.1, 1 + 0.1],
        tickmode="array",
        row=4,
        col=1,
        secondary_y=False,
    )
    fig.update_xaxes(
        title="Position", row=4, col=1
    )  # Add the x-axis title to the bottom plot

    # fig.update_layout(paper_bgcolor="#000000")
    fig.update_layout(height=900, title_text=chromname, showlegend=False)
    # plot_bgcolor="rgba(246, 237, 202, 0.5)",

    # Remove unnecessary plot elements
    fig.update_layout(
        modebar_remove=[
            "select2d",
            "lasso2d",
            "select",
            "autoScale2d",
            "zoom",
            "toImage",
        ]
    )

    # Write a png version of the plot
    fig.write_image(
        str(outdir.absolute() / (chromname + ".png")),
        format="png",
        height=900,
        width=1600,
    )

    # Write a html version of the plot
    return fig.to_html(
        include_plotlyjs=True if offline_plots else "cdn", full_html=False
    )
