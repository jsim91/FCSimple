library(shinythemes)
library(flowCore)
library(flowWorkspace)

# to do: store final transform values, save out transform parameters;
# allow user to update these values every time the button is pressed; finally, apply transform to data set
# and produce an output with the transform values used, allow an option to parse this output for future
# transform setting on files using identical panels, ensures consistency and streamlines workflow consider
# a json format or a csv file with columns channel names, and rows parameters, perhaps a row for alternate
# names, such as fluorochrome + marker or metal + marker to allow for more robust matching of future panels

# Define UI ----
fluidPage(theme = shinytheme("yeti"),
                navbarPage(
                  "FCSimple Transformer",
                  tabPanel("Transform",
                           sidebarPanel(
                             uiOutput("channel"),
                             uiOutput("transform_type"),
                             uiOutput("hyperparameters"),
                             hr(),
                             actionButton(inputId = "apply_transform", label = "apply transform to channel"),
                             hr(),
                             actionButton(inputId = "apply_transform_all", label = "apply transform to all"),
                             hr(),
                             uiOutput("channel_x"),
                             uiOutput("channel_y"),
                             actionButton(inputId = "apply_plot", label = "update 2d plot")
                           ),
                           mainPanel(
                             plotOutput("densPlot"),
                             textOutput("dens_hyperlog_error"),
                             textOutput("dens_asinh_error"),
                             hr(),
                             plotOutput("render_2d_plot")
                           )),
                  tabPanel("Selected Transforms",
                           tableOutput("parameter_df")),
                  tabPanel("UMAP Preview",
                           sidebarPanel(
                             uiOutput("init_input"),
                             uiOutput("min_dist_input"),
                             uiOutput("spread_input"),
                             uiOutput("neighbors_input"),
                             actionButton(inputId = "umap_button", label = "calculate sample UMAP"),
                           ),
                           mainPanel(
                             plotOutput("sample_umap")
                           )),
                  tabPanel("Finalize",
                           mainPanel(
                             column(12, align = "center",
                               uiOutput("warning_message"),
                               hr(),
                               actionButton(inputId = "finalize_transform", label = "Yes, finalize transformation choices")
                             )
                           ))
                )
)
