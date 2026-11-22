project <- path.expand("~/Documents/st773_epidemiology_96_cRMURn")

if (!dir.exists(project)) {
  archive <- paste0(project, ".tar.gz")
  stopifnot(file.exists(archive))
  dir.create(project)
  stopifnot(untar(archive, exdir = project) == 0)
}

timeline_file <- list.files(
  project, "^st773_collection_timeline\\.tsv$",
  recursive = TRUE, full.names = TRUE
)
stopifnot(length(timeline_file) == 1L)

setwd(dirname(dirname(timeline_file)))
timeline <- read.delim(timeline_file, colClasses = "character")
print(names(timeline))
print(head(timeline, 2))

additions <- read.delim(
  text = "accession|biosample|bioproject|published_isolate|collection_date|reported_week|collection_location|date_source
GCA_025891215.1|SAMN31309206||D2|2022-03 to 2022-06||Germany: Facility 5|LeBreton et al. supplement
GCA_027130495.1|SAMN31807497||RIVM_C046766|2022-07-15||Netherlands|LeBreton et al. supplement
GCA_027132365.1|SAMN31807506||RIVM_C047413|2022-08-05||Netherlands|LeBreton et al. supplement
GCA_027131475.1|SAMN31807510||RIVM_C047473|2022-08-18||Netherlands|LeBreton et al. supplement
GCA_027131535.1|SAMN31807521||RIVM_C048075|2022-08-24||Netherlands|LeBreton et al. supplement
GCA_027908815.1|SAMEA112251230||S.161.22.Pa|2022-08-22||Germany: Facility 3|LeBreton et al. supplement
GCA_027908735.1|SAMEA112251232||S.207.22.Pa|2022-10-10||Germany: Facility 3|LeBreton et al. supplement
GCA_029623785.1|SAMN33993894||110817|2022-11-26||Germany: Facility 2|LeBreton et al. supplement
GCA_062088240.1|SAMN60549723||25|2025-04-14||Ukraine|NCBI BioSample
GCA_062088236.1|SAMN60549726||32|2025-06-04||Ukraine|NCBI BioSample
ERZ26343803|SAMEA113542066|PRJEB88251|9953|2022-01-17||Greece: Attica H1|Attica study metadata
ERZ26343808|SAMEA113542067|PRJEB88251|10071|2022-12-12||Greece: Attica H2|Attica study metadata
ERZ26343812|SAMEA113542068|PRJEB88251|9912|2022-10-06||Greece: Attica H2|Attica study metadata
ERZ26343816|SAMEA113542069|PRJEB88251|9958|2022-05-17||Greece: Attica H1|Attica study metadata
ERZ26343821|SAMEA113542070|PRJEB88251|9964|2022-07-23||Greece: Attica H1|Attica study metadata
ERR11021105|SAMEA112785803|PRJEB60743|AMA004593|||Denmark|
ERR11021118|SAMEA112785801|PRJEB60743|AMA004509|||Denmark|
ERR11021131|SAMEA112785802|PRJEB60743|AMA004378|||Denmark|
ERR11021143|SAMEA112785800|PRJEB60743|AMA004584|||Denmark|
GCA_003725635.1|||PS1|||Hungary|",
  sep = "|", colClasses = "character", quote = ""
)

timeline96 <- rbind(timeline, additions)

snp <- read.delim(
  "results/st773_parsnp_v12_96/st773_snp_distances.tsv",
  row.names = 1, check.names = FALSE
)
matrix_ids <- sub("\\.fna(\\.ref)?$", "", basename(rownames(snp)))
stopifnot(
  nrow(timeline96) == 96L,
  length(matrix_ids) == 96L,
  anyDuplicated(timeline96$accession) == 0L,
  setequal(timeline96$accession, matrix_ids)
)

timeline96 <- timeline96[match(matrix_ids, timeline96$accession), ]
timeline96$country <- trimws(sub(":.*$", "", timeline96$collection_location))
timeline96$country[timeline96$country %in% "Czech Republic"] <- "Czechia"

out <- "metadata/st773_figure1_timeline_96.tsv"
write.table(timeline96, out, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

has_date <- !is.na(timeline96$collection_date) &
  nzchar(trimws(timeline96$collection_date))
print(table(
  Country = timeline96$country,
  Date = ifelse(has_date, "Recorded", "Missing"),
  useNA = "ifany"
))
cat("Saved", nrow(timeline96), "assemblies to", out, "\n")

#Figure1Α
year_ticks <- as.Date(paste0(
  c(2009, 2012, 2015, 2018, 2020, 2022, 2024, 2025), "-01-01"
))

fig1a_polished <- fig1a +
  geom_segment(
    data = anchors,
    aes(x = label_x, xend = label_x, y = y + 0.09, yend = label_y - 0.04,
        colour = country),
    inherit.aes = FALSE, linewidth = 0.25, alpha = 0.6
  ) +
  scale_x_date(breaks = year_ticks, date_labels = "%Y") +
  labs(
    title = "A   Sampling chronology",
    subtitle = "96 ST773 assemblies from nine European countries",
    x = "Collection year",
    caption = paste(
      "Points show exact collection dates. Bars show reported date ranges.",
      "Dates unavailable for seven assemblies (Denmark 4, Hungary 1, Spain 2).",
      sep = "\n"
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(colour = "#758392", linewidth = 0.35),
    axis.text = element_text(colour = "#2E3E4B"),
    plot.title = element_text(face = "bold", size = 13, colour = "#182C3B"),
    plot.subtitle = element_text(colour = "#526373", margin = margin(b = 8)),
    plot.caption = element_text(hjust = 0, size = 8, colour = "#526373"),
    plot.caption.position = "plot",
    plot.margin = margin(9, 12, 8, 9)
  )

ggsave("figures/figure1a_sampling_timeline_polished.svg", fig1a_polished,
       device = svglite::svglite, width = 180, height = 108,
       units = "mm", bg = "white")
print(fig1a_polished)

print(subset(
  dated, country == "Norway" & format(start, "%Y") == "2009",
  select = c(accession, biosample, published_isolate, collection_date, date_source)
), row.names = FALSE)

# Figure 1B

packages <- c("mapsf", "sf", "svglite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

# Basemap from QGIS, derived from Natural Earth
basemap_file <- file.path(tempdir(), "qgis_world_map.gpkg")
if (!file.exists(basemap_file)) {
  download.file(
    "https://raw.githubusercontent.com/qgis/QGIS/master/resources/data/world_map.gpkg",
    basemap_file, mode = "wb"
  )
}

world <- sf::st_read(
  basemap_file, layer = "countries", quiet = TRUE,
  wkt_filter = "POLYGON ((-15 32, 52 32, 52 64, -15 64, -15 32))"
)
world <- sf::st_make_valid(sf::st_set_crs(world, NA))
europe <- suppressWarnings(sf::st_crop(
  world, xmin = -15, xmax = 52, ymin = 32, ymax = 64
))
europe <- sf::st_transform(sf::st_set_crs(europe, 4326), 3035)
involved <- europe$iso_a2 %in% c("ES", "CZ", "NL", "DE", "PL", "TR", "NO", "DK")
ukraine <- europe$iso_a2 %in% "UA"

sites <- data.frame(
  name = c("Barcelona", "Pilsen", "Gelderland", "Landstuhl",
           "Ukraine", "Kyiv", "Poland", "Türkiye"),
  lon = c(2.17, 13.38, 5.94, 7.57, 29.15, 30.52, 19.40, 34.30),
  lat = c(41.39, 49.75, 52.08, 49.42, 48.20, 50.45, 52.05, 39.30),
  dx = c(0, 0, -60000, -100000, 90000, 70000, 0, 0),
  dy = c(-145000, -135000, 150000, -65000, -105000, 80000, 135000, -140000),
  align = c(.5, .5, 1, 1, 0, 0, .5, .5)
)
xy <- sf::st_coordinates(sf::st_transform(
  sf::st_as_sf(sites, coords = c("lon", "lat"), crs = 4326), 3035
))
rownames(xy) <- sites$name

links <- data.frame(
  from = c("Ukraine", "Türkiye", "Ukraine", "Ukraine", "Poland", "Kyiv"),
  to = c("Barcelona", "Ukraine", "Pilsen", "Poland", "Gelderland", "Landstuhl"),
  case = c(1, 2, 2, 3, 3, 4),
  bend = c(-.18, -.11, -.10, -.12, .22, .12),
  transfer = c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE),
  arrowhead = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
)

colours <- c("#A33E68", "#2166AC", "#B05C25", "#39724C")
destinations <- c(
  "Barcelona, Spain", "Pilsen, Czechia",
  "Gelderland, Netherlands", "Landstuhl, Germany"
)
case_notes <- c(
  "Surgery in Ukraine\nAdmitted in February\nRectal isolate, 14 April",
  "Transplant in Türkiye (2021)\nFollow-up care in Ukraine\nMay repatriation and detection",
  "Ukrainian soldier\nTransferred from a Polish\nmilitary hospital\nRectal screening swab, week 34",
  "Care near Dnipro, then in Kyiv\n26 November specimen"
)

# Schematic country markers for the Norwegian and Danish surveillance links
cohorts <- data.frame(
  name = c("Norway", "Denmark"),
  lon = c(9.05, 9.32),
  lat = c(60.90, 56.13),
  bend = c(.13, .05),
  note = c(
    "Prior hospitalisation in Ukraine, 2022",
    "Ukraine-linked surveillance cohort, 2022–23"
  )
)
cohort_xy <- sf::st_coordinates(sf::st_transform(
  sf::st_as_sf(cohorts, coords = c("lon", "lat"), crs = 4326), 3035
))

curved_path <- function(start, end, bend) {
  delta <- end - start
  distance <- sqrt(sum(delta^2))
  control <- (start + end) / 2 + bend * c(delta[2], -delta[1])
  t <- seq(30000 / distance, 1 - 75000 / distance, length.out = 100)
  outer((1 - t)^2, start) + outer(2 * t * (1 - t), control) + outer(t^2, end)
}

draw_healthcare_map <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  ink <- "#152737"
  muted <- "#4A6070"
  cohort_ink <- "#536C7C"
  
  par(fig = c(0, 1, 0, 1), mar = rep(0, 4), oma = rep(0, 4),
      xaxs = "i", yaxs = "i", family = "sans", ps = 10, cex = 1)
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  
  text(.028, .944, "B", adj = 0, font = 2, cex = 1.3, col = ink)
  text(.067, .944, "Cross-border healthcare histories",
       adj = 0, font = 2, cex = 1.3, col = ink)
  text(.067, .891, "NDM-1-producing ST773 detections linked to Ukraine, 2022–23",
       adj = 0, cex = .88, col = muted)
  
  key_y <- c(.803, .676, .549, .402)
  points(rep(.685, 4), key_y, pch = 21, bg = colours, col = NA, cex = 1.35)
  text(rep(.685, 4), key_y, 1:4, cex = .72, font = 2, col = "white")
  text(rep(.710, 4), key_y, destinations, adj = 0, font = 2, cex = .86, col = ink)
  
  for (i in seq_along(case_notes)) {
    note_lines <- strsplit(case_notes[i], "\n", fixed = TRUE)[[1]]
    note_y <- key_y[i] - .033 - (seq_along(note_lines) - 1) * .024
    text(rep(.710, length(note_lines)), note_y, note_lines,
         adj = 0, cex = .74, col = muted)
  }
  
  separators <- c(.701, .574, .426, .306)
  segments(.676, separators, .969, separators, col = "#D7E0E4")
  
  cohort_y <- c(.275, .193)
  points(rep(.685, 2), cohort_y, pch = 22, bg = "white", col = cohort_ink, cex = .86)
  text(rep(.710, 2), cohort_y, cohorts$name, adj = 0, font = 2, cex = .86, col = ink)
  text(rep(.710, 2), cohort_y - .029, cohorts$note, adj = 0, cex = .73, col = muted)
  
  arrows(.035, .094, .077, .094, length = .07, lwd = 1.5, col = muted)
  text(.085, .094, "Documented patient transfer", adj = 0, cex = .75, col = ink)
  segments(.327, .094, .369, .094, lty = 2, lwd = 1.5, col = muted)
  text(.377, .094, "Reported care or evacuation history", adj = 0, cex = .75, col = ink)
  segments(.681, .094, .723, .094, lty = 3, lwd = 1.4, col = cohort_ink)
  text(.731, .094, "Country-level surveillance link", adj = 0, cex = .75, col = ink)
  text(.035, .047, "Paths and country-level endpoints are schematic. Basemap from Natural Earth.",
       adj = 0, cex = .71, col = muted)
  
  par(fig = c(.027, .661, .145, .839), mar = rep(0, 4), new = TRUE)
  plot.new()
  plot.window(c(2.80, 6.69) * 1e6, c(1.66, 4.43) * 1e6, asp = 1)
  limits <- par("usr")
  rect(limits[1], limits[3], limits[2], limits[4], col = "#EEF4F6", border = NA)
  
  mapsf::mf_map(europe, col = "#F5F6F4", border = "#C3CCD2", lwd = .55, add = TRUE)
  mapsf::mf_map(europe[involved, ], col = "#DEE9EE",
                border = "#C3CCD2", lwd = .55, add = TRUE)
  mapsf::mf_map(europe[ukraine, ], col = "#CCDDE8",
                border = "#C3CCD2", lwd = .55, add = TRUE)
  
  for (i in seq_len(nrow(cohorts))) {
    path <- curved_path(xy["Ukraine", ], cohort_xy[i, ], cohorts$bend[i])
    lines(path, col = cohort_ink, lwd = 1.4, lty = 3)
  }
  
  for (i in seq_len(nrow(links))) {
    path <- curved_path(xy[links$from[i], ], xy[links$to[i], ], links$bend[i])
    colour <- colours[links$case[i]]
    line_type <- if (links$transfer[i]) 1 else 2
    
    lines(path, col = "white", lwd = 3.5, lty = line_type)
    lines(path, col = colour, lwd = 1.9, lty = line_type)
    
    if (links$arrowhead[i]) {
      arrows(path[93, 1], path[93, 2], path[100, 1], path[100, 2],
             length = .08, angle = 25, lwd = 1.9, col = colour)
    }
  }
  
  points(xy[5:8, ], pch = 21, bg = "white", col = muted, cex = .65)
  points(xy[1:4, ], pch = 21, bg = colours, col = "white", cex = 1.55)
  text(xy[1:4, ], labels = 1:4, col = "white", font = 2, cex = .73)
  
  for (i in seq_len(nrow(sites))) {
    text(xy[i, 1] + sites$dx[i], xy[i, 2] + sites$dy[i], sites$name[i],
         adj = c(sites$align[i], .5), cex = .74, col = ink)
  }
  
  points(cohort_xy, pch = 22, bg = "white", col = cohort_ink, cex = .86)
  text(cohort_xy[, 1] + 65000, cohort_xy[, 2] + 20000, cohorts$name,
       adj = 0, cex = .74, col = ink)
  box(col = "#CAD4DA", lwd = .6)
}

dir.create("figures", showWarnings = FALSE)
stem <- "figures/figure1b_st773_healthcare_map"

svglite::svglite(paste0(stem, ".svg"), width = 180 / 25.4, height = 108 / 25.4)
invisible(tryCatch(draw_healthcare_map(), finally = dev.off()))

png(paste0(stem, ".png"), width = 180, height = 108, units = "mm", res = 600)
invisible(tryCatch(draw_healthcare_map(), finally = dev.off()))

draw_healthcare_map()

#Figure2A
packages <- c("ape", "svglite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)

tree_file <- "results/st773_parsnp_v12_96/iqtree2/st773_96.treefile"
metadata_file <- "metadata/st773_figure1_timeline_96.tsv"
stopifnot(file.exists(tree_file), file.exists(metadata_file))

tree <- ape::ladderize(ape::read.tree(tree_file))
meta <- read.delim(metadata_file, colClasses = "character")

madrid <- meta$accession %in% c("GCA_036321415.2", "GCA_036321435.2")
stopifnot(sum(madrid) == 2L)

missing_year <- madrid & (
  is.na(meta$collection_date) | trimws(meta$collection_date) == ""
)

if (any(missing_year)) {
  meta$collection_date[missing_year] <- "2022"
  meta$date_source[missing_year] <- "Pitart et al. (2025), Figure 1"
  write.table(meta, metadata_file, sep = "\t", row.names = FALSE,
              quote = FALSE, na = "")
}

accession <- function(x) {
  pattern <- "((?:GCA|GCF)_[0-9]+\\.[0-9]+|(?:ERR|SRR|DRR|ERZ)[0-9]+)"
  sub(paste0("^.*?", pattern, ".*$"), "\\1", x, perl = TRUE)
}

tip_id <- accession(tree$tip.label)
meta_id <- accession(meta$accession)
stopifnot(length(tip_id) == 96L, !anyDuplicated(tip_id), !anyDuplicated(meta_id))

matched <- match(tip_id, meta_id)
if (anyNA(matched)) {
  stop("Tree tips missing from metadata: ", paste(tip_id[is.na(matched)], collapse = ", "))
}
meta <- meta[matched, , drop = FALSE]

palette <- c(
  Czechia = "#2166AC", Denmark = "#687384", Germany = "#39724C",
  Greece = "#927521", Hungary = "#8B644E", Netherlands = "#B05C25",
  Norway = "#137F89", Spain = "#A33E68", Ukraine = "#60489A"
)

location <- meta$collection_location
if ("country" %in% names(meta)) {
  known <- !is.na(meta$country) & nzchar(trimws(meta$country))
  location[known] <- meta$country[known]
}
location <- sub("Czech Republic", "Czechia", location, ignore.case = TRUE)

country <- rep(NA_character_, nrow(meta))
for (name in names(palette)) {
  country[grepl(name, location, ignore.case = TRUE)] <- name
}
if (anyNA(country)) {
  stop("Check country for: ", paste(tip_id[is.na(country)], collapse = ", "))
}

year_pattern <- "(?:19|20)[0-9]{2}"
year <- ifelse(
  grepl(year_pattern, meta$collection_date),
  sub(paste0("^.*?(", year_pattern, ").*$"), "\\1", meta$collection_date, perl = TRUE),
  "\u2013"
)

anchors <- c(
  GCA_043266445.1 = "Ukraine, Nov 2019",
  GCA_062088240.1 = "Ukraine, Apr 2025",
  ERR13374403 = "Barcelona, Apr 2022",
  GCA_036872625.1 = "Pilsen, May 2022"
)
anchor_rows <- match(names(anchors), tip_id)
stopifnot(!anyNA(anchor_rows))

support <- strsplit(as.character(tree$node.label), "/", fixed = TRUE)
supported_nodes <- length(tip_id) + which(vapply(support, function(pair) {
  if (length(pair) != 2L) return(FALSE)
  values <- suppressWarnings(as.numeric(pair))
  all(!is.na(values) & values >= c(80, 95))
}, logical(1)))

draw_tree <- function() {
  old_par <- par(mar = c(0.7, 0.8, 0.7, 0.8), family = "sans", ps = 10)
  on.exit(par(old_par))
  
  n <- length(tip_id)
  span <- max(ape::node.depth(tree, method = 2)) - 1
  
  ape::plot.phylo(
    tree, use.edge.length = FALSE, node.depth = 2, show.tip.label = FALSE,
    edge.color = "#344555", edge.width = 0.8,
    x.lim = c(-0.02, 1.60) * span, y.lim = c(-10, n + 8)
  )
  layout <- get("last_plot.phylo", envir = getFromNamespace(".PlotPhyloEnv", "ape"))
  y <- layout$yy[seq_len(n)]
  
  points(layout$xx[supported_nodes], layout$yy[supported_nodes],
         pch = 16, cex = 0.42, col = "#344555")
  
  rect(1.035 * span, y - 0.41, 1.085 * span, y + 0.41,
       col = unname(palette[country]), border = NA)
  text(1.155 * span, y, year, cex = 0.57, col = "#344555")
  
  text(0, n + 5.5, "A   Core-SNP phylogeny", adj = 0, font = 2, cex = 1.18)
  text(1.060 * span, n + 2.2, "Country", font = 2, cex = 0.70)
  text(1.155 * span, n + 2.2, "Year", font = 2, cex = 0.70)
  
  ordered <- order(y[anchor_rows])
  rows <- anchor_rows[ordered]
  label_y <- y[rows]
  for (i in seq.int(2L, length(label_y))) {
    label_y[i] <- max(label_y[i], label_y[i - 1L] + 2.3)
  }
  label_y <- label_y - max(0, max(label_y) - n)
  
  segments(1.20 * span, y[rows], 1.245 * span, label_y, col = "#9CA9B2")
  text(1.255 * span, label_y, unname(anchors[ordered]), adj = 0,
       font = 2, cex = 0.74, col = unname(palette[country[rows]]))
  
  legend(0, -0.4, legend = names(palette), fill = unname(palette),
         border = NA, bty = "n", ncol = 5, cex = 0.73, x.intersp = 0.6)
  points(0, -7.4, pch = 16, cex = 0.52, col = "#344555")
  text(0.025 * span, -7.4,
       "SH-aLRT \u2265 80% and UFBoot \u2265 95%. Branch lengths are schematic.",
       adj = 0, cex = 0.67, col = "#344555")
}

dir.create("figures", showWarnings = FALSE)
output <- "figures/figure2a_st773_core_snp_tree"

svglite::svglite(paste0(output, ".svg"), width = 180 / 25.4, height = 210 / 25.4)
invisible(tryCatch(draw_tree(), finally = dev.off()))

png(paste0(output, ".png"), width = 180, height = 210, units = "mm", res = 600)
invisible(tryCatch(draw_tree(), finally = dev.off()))

draw_tree()
cat("Matched", length(tip_id), "tree tips and saved SVG and 600 dpi PNG.\n")


#Figure2B
packages <- c("beeswarm", "svglite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)

matrix_file <- "results/st773_parsnp_v12_96/st773_snp_distances.tsv"
metadata_file <- "metadata/st773_figure1_timeline_96.tsv"

accession <- function(x) {
  pattern <- "((?:GCA|GCF)_[0-9]+\\.[0-9]+|(?:ERR|SRR|DRR|ERZ)[0-9]+)"
  sub(paste0("^.*?", pattern, ".*$"), "\\1", x, perl = TRUE)
}

distances <- as.matrix(read.delim(matrix_file, row.names = 1, check.names = FALSE))
rownames(distances) <- accession(rownames(distances))
colnames(distances) <- accession(colnames(distances))
stopifnot(
  identical(dim(distances), c(96L, 96L)),
  is.numeric(distances), all(is.finite(distances)),
  !anyDuplicated(rownames(distances)), !anyDuplicated(colnames(distances)),
  setequal(rownames(distances), colnames(distances))
)
distances <- distances[, rownames(distances), drop = FALSE]
stopifnot(isSymmetric(distances), all(diag(distances) == 0))

focal <- "GCA_062088240.1"
stopifnot(focal %in% rownames(distances))
ids <- setdiff(colnames(distances), focal)
snp <- distances[focal, ids]

meta <- read.delim(metadata_file, colClasses = "character")
meta_ids <- accession(meta$accession)
matched <- match(ids, meta_ids)
stopifnot(!anyDuplicated(meta_ids), !anyNA(matched))
meta <- meta[matched, , drop = FALSE]

colours <- c(
  Czechia = "#2166AC", Denmark = "#687384", Germany = "#39724C",
  Greece = "#927521", Hungary = "#8B644E", Netherlands = "#B05C25",
  Norway = "#137F89", Spain = "#A33E68", Ukraine = "#60489A"
)
location <- meta$collection_location
if ("country" %in% names(meta)) {
  known <- !is.na(meta$country) & nzchar(trimws(meta$country))
  location[known] <- meta$country[known]
}
location <- sub("Czech Republic", "Czechia", location, ignore.case = TRUE)

country <- rep(NA_character_, length(ids))
for (name in names(colours)) {
  country[grepl(name, location, ignore.case = TRUE)] <- name
}
if (anyNA(country)) {
  stop("Check country for: ", paste(ids[is.na(country)], collapse = ", "))
}

selected <- c(
  GCA_043266445.1 = "Nov 2019",
  GCA_036872625.1 = "Pilsen, May 2022",
  GCA_043265905.1 = "Aug 2022",
  GCA_041931795.1 = "Dec 2022",
  GCA_046498825.1 = "PAN-58, 2022",
  ERR13374403 = "Barcelona, Apr 2022"
)
stopifnot(all(names(selected) %in% ids))

top_to_bottom <- c("Ukraine", setdiff(names(colours), "Ukraine"))
groups <- split(ids, factor(country, levels = rev(top_to_bottom)))
plot_ids <- unlist(groups, use.names = FALSE)
swarm_values <- lapply(groups, function(group_ids) log10(1 + snp[group_ids]))

draw_distances <- function() {
  old_par <- par(mar = c(3.2, 7.4, 2.7, 1.2), ps = 10, family = "sans")
  on.exit(par(old_par))
  
  ticks <- c(0, 5, 10, 20, 50, 100, 250, 500, 1000, 2500, 5000)
  if (min(snp) >= 5) ticks <- ticks[ticks > 0]
  tick_x <- log10(1 + ticks)
  x_limits <- range(tick_x) + c(-0.05, 0.05)
  
  plot.new()
  plot.window(xlim = x_limits, ylim = c(0, 1), xaxs = "i", yaxs = "i")
  
  # Measure each swarm before assigning its row height.
  n_groups <- length(groups)
  group_index <- rep(seq_len(n_groups), lengths(groups))
  swarm <- beeswarm::beeswarm(
    swarm_values, horizontal = TRUE, add = TRUE, do.plot = FALSE,
    at = seq_len(n_groups), cex = 0.65
  )
  offsets <- swarm$x - group_index
  swarm_height <- 2 * as.numeric(tapply(abs(offsets), group_index, max))
  shrink <- min(1, 0.65 / max(sum(swarm_height), 0.001))
  row_height <- shrink * swarm_height
  row_height <- row_height + (1 - sum(row_height)) / n_groups
  edges <- c(0, cumsum(row_height))
  centres <- head(edges, -1) + row_height / 2
  point_y <- centres[group_index] + shrink * offsets
  
  shaded <- seq.int(1L, n_groups, by = 2L)
  rect(x_limits[1], edges[shaded], x_limits[2], edges[shaded + 1L],
       col = "#F5F7F9", border = NA)
  segments(tick_x, 0, tick_x, 1, col = "#E3E8EC", lwd = 0.7)
  
  point_colours <- unname(colours[names(groups)[group_index]])
  marked <- plot_ids %in% names(selected)
  points(swarm$y[!marked], point_y[!marked], pch = 16,
         col = point_colours[!marked], cex = 0.65 * shrink)
  points(swarm$y[marked], point_y[marked], pch = 23, bg = "white",
         col = point_colours[marked], cex = 0.94 * shrink, lwd = 1.3)
  
  axis(2, at = centres, labels = sprintf("%s (%d)", names(groups), lengths(groups)),
       las = 1, tick = FALSE, cex.axis = 0.89, col.axis = "#223646")
  axis(1, at = tick_x, labels = format(ticks, big.mark = ",", trim = TRUE),
       cex.axis = 0.82, tcl = -0.2, mgp = c(2, 0.55, 0),
       col = "#8B99A5", col.axis = "#223646")
  
  mtext("B   SNP distances to the April 2025 Ukrainian isolate",
        side = 3, line = 1.0, adj = 0, font = 2, cex = 1.13)
  mtext("Core SNP differences", side = 1, line = 1.95, cex = 0.92)
}
dir.create("figures", showWarnings = FALSE)
output <- "figures/figure2b_st773_april2025_snp_distances"

svglite::svglite(paste0(output, ".svg"), width = 180 / 25.4, height = 155 / 25.4)
invisible(tryCatch(draw_distances(), finally = dev.off()))

png(paste0(output, ".png"), width = 180, height = 155, units = "mm", res = 600)
invisible(tryCatch(draw_distances(), finally = dev.off()))

draw_distances()

print(data.frame(
  comparator = unname(selected),
  accession = names(selected),
  core_SNPs = unname(snp[names(selected)])
), row.names = FALSE)
