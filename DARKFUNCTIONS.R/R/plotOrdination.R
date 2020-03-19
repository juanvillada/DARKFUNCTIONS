#' @export

plotOrdination <- function(DATA, METADATA, MAPPING) {
  attach(METADATA, warn.conflicts = F)

  if (MAPPING == "Layer") {
    COLORS <- c("#b7d8ff", "#98c699")
  }
  if (MAPPING == "Vegetation") {
    COLORS <- c("#dfc3f8", "#beefc1", "#eca6c1", "#61b7d9", "#f9b99f")
  }
  if (MAPPING == "Ecosystem") {
    COLORS <- c("#dfc3f8", "#eca6c1", "#f9b99f")
  }

  plot(DATA, display = "sites", type = "n")

  if (MAPPING == "Layer") {
    points(DATA, display = "sites", pch = c(16, 17)[Layer], col = COLORS[Layer])
    ordiellipse(DATA, groups = Layer, kind = "sd", draw = "polygon", col = COLORS)
    legend("bottomleft", legend = levels(metadata$Layer), bty = "n", col = COLORS, pch = c(16, 17))
  }
  if (MAPPING == "Vegetation") {
    points(DATA, display = "sites", pch = c(16, 17)[Layer], col = COLORS[Vegetation])
    ordiellipse(DATA, groups = Vegetation, kind = "sd", draw = "polygon", col = COLORS)
    legend("bottomleft", legend = levels(metadata$Layer), bty = "n", pch = c(16, 17))
    legend("bottomright", legend = levels(metadata$Vegetation), bty = "n", col = COLORS, pch =15)
  }
  if (MAPPING == "Ecosystem") {
    points(DATA, display = "sites", pch = c(16, 17)[Layer], col = COLORS[Ecosystem])
    ordiellipse(DATA, groups = Ecosystem, kind = "sd", draw = "polygon", col = COLORS)
    legend("bottomleft", legend = levels(metadata$Layer), bty = "n", pch = c(16, 17))
    legend("bottomright", legend = levels(metadata$Ecosystem), bty = "n", col = COLORS, pch =15)
  }

  detach(METADATA)
}
