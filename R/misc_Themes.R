# Helper functions --------------------------------------------------------

theme_Job <- ggplot2::theme(
    text = ggplot2::element_text(size = 8, family='Nimbus Sans', face = 'bold'),
    axis.text = ggtext::element_markdown(family='Nimbus Sans', face = 'bold'),
    axis.text.x = ggtext::element_markdown(family='Nimbus Sans', face = 'bold'),
    axis.title.x = ggtext::element_textbox_simple(family='Nimbus Sans', face = 'bold', width = NULL, halign = .5),
    axis.title.y = ggtext::element_textbox_simple(family='Nimbus Sans', face = 'bold', orientation = 'left-rotated', width = NULL, size = 8, halign = .5),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(colour = 'grey75', linetype = '12'),
    panel.grid.minor.y = ggplot2::element_line(colour = 'grey75', linetype = '12'),
    panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
    panel.border = ggplot2::element_rect(fill = NA, colour = NA),
    legend.position = 'bottom',
    legend.text = ggtext::element_markdown(family='Nimbus Sans', face = 'bold'),
    legend.title = ggtext::element_markdown(family='Nimbus Sans', face = 'bold'),
)

theme_Job2 <- ggplot2::theme(
    text = ggplot2::element_text(size = 8, family='Nimbus Sans', face = 'bold'),
    axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
    axis.title.x = ggtext::element_textbox_simple(family='Nimbus Sans', face = 'bold', width = NULL, halign = .5),
    axis.title.y = ggtext::element_textbox_simple(family='Nimbus Sans', face = 'bold', orientation = 'left-rotated', width = NULL, size = 8, halign = .5),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
    panel.border = ggplot2::element_rect(fill = NA, colour = NA)
)
