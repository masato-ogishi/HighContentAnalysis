singleCellCropping <- function(cellID, dataDir="./Image_Train/", debriMinSizeThr=100){
  
  ## Check directory
  dir.create(file.path(dataDir, "Original"), recursive=T, showWarnings=F)
  dir.create(file.path(dataDir, "Segmented"), recursive=T, showWarnings=F)
  dir.create(file.path(dataDir, "Cropped"), recursive=T, showWarnings=F)
  
  ## Load images
  originalHeader <- file.path(dataDir, "Original", cellID)
  r <- EBImage::readImage(paste0(originalHeader, "_red.png"))
  g <- EBImage::readImage(paste0(originalHeader, "_green.png"))
  b <- EBImage::readImage(paste0(originalHeader, "_blue.png"))
  y <- EBImage::readImage(paste0(originalHeader, "_yellow.png"))
  img <- EBImage::rgbImage(r, g, b)
  halfSize <- dim(img)[1]/2
  
  ## Nucleus segmentation using the Blue channel
  nmask <- EBImage::thresh(b, w=round(dim(img)[1]/20), h=round(dim(img)[2]/20))
  nseg <- EBImage::bwlabel(nmask)
  nf <- EBImage::computeFeatures.shape(nseg)
  nr <- which(nf[,"s.area"] < debriMinSizeThr)
  nseg <- EBImage::rmObjects(nseg, nr)
  nseg <- EBImage::closing(nseg, EBImage::makeBrush(15, shape='diamond'))
  nseg <- EBImage::fillHull(nseg)
  nseg <- EBImage::watershed(EBImage::distmap(nseg), tolerance=1)
  if(identical(max(nseg), 0L)) return(NULL)
  
  ## Membrane segmentation using the Red and Green channels
  rg.sqrt <- sqrt(r^2+g^2)
  cmask <- EBImage::filter2(rg.sqrt, EBImage::makeBrush(7, shape='diamond'))
  cmask <- EBImage::thresh(cmask, w=round(dim(img)[1]/80), h=round(dim(img)[2]/80))
  cseg <- EBImage::bwlabel(cmask)
  cf <- EBImage::computeFeatures.shape(cseg)
  cr <- which(cf[,"s.area"] < debriMinSizeThr)
  cseg <- EBImage::rmObjects(cseg, cr)
  cseg <- EBImage::dilate(cseg, EBImage::makeBrush(15, shape='diamond'))
  cseg <- EBImage::fillHull(cseg)
  cseg <- (cseg + nseg)
  cseg[cseg>=1] <- 1L
  rgb.sqrt <- sqrt(r^2+g^2+b^2)
  cseg <- EBImage::propagate(rgb.sqrt, seeds=nseg, mask=cseg)
  cseg <- EBImage::fillHull(cseg)
  if(identical(max(cseg), 0L)) return(NULL)
  
  ### Cropping individual cells [+ clipping]
  #outHeader1 <- file.path(dataDir, "Cropped", cellID)
  #cropSize <- mean(dim(img)[1:2])/4
  #cropSize.half <- floor(cropSize/2)
  #cf <- data.table::as.data.table(EBImage::computeFeatures.moment(cseg))
  #cf[,CellID:=.I]
  #cf[,m.cx:=floor(m.cx)][,m.cy:=floor(m.cy)]
  #cf[,X1:=m.cx-cropSize.half][,X2:=m.cx+cropSize.half][,Y1:=m.cy-cropSize.half][,Y2:=m.cy+cropSize.half]
  #cf <- cf[X1>=1,][X2<=dim(img)[1],][Y1>=1,][Y2<=dim(img)[2],]
  #lapply(cf$CellID, function(i){
  #  cseg.crop <- cseg==i
  #  cropRange <- unlist(cf[CellID==i, .(X1, X2, Y1, Y2)])
  #  r.crop <- EBImage::resize((r*cseg.crop)[cropRange["X1"]:cropRange["X2"],cropRange["Y1"]:cropRange["Y2"]], w=100, h=100)
  #  g.crop <- EBImage::resize((g*cseg.crop)[cropRange["X1"]:cropRange["X2"],cropRange["Y1"]:cropRange["Y2"]], w=100, h=100)
  #  b.crop <- EBImage::resize((b*cseg.crop)[cropRange["X1"]:cropRange["X2"],cropRange["Y1"]:cropRange["Y2"]], w=100, h=100)
  #  y.crop <- EBImage::resize((y*cseg.crop)[cropRange["X1"]:cropRange["X2"],cropRange["Y1"]:cropRange["Y2"]], w=100, h=100)
  #  EBImage::writeImage(r.crop, paste0(outHeader1, "_red_crop_", i, ".png"))
  #  EBImage::writeImage(g.crop, paste0(outHeader1, "_green_crop_", i, ".png"))
  #  EBImage::writeImage(b.crop, paste0(outHeader1, "_blue_crop_", i, ".png"))
  #  EBImage::writeImage(y.crop, paste0(outHeader1, "_yellow_crop_", i, ".png"))
  #})
  
  ## Cropping individual cells [clipping(-)]
  outHeader1 <- file.path(dataDir, "Cropped", cellID)
  cf <- data.table::as.data.table(EBImage::computeFeatures.moment(cseg))
  cf[,CellID:=.I]
  lapply(cf$CellID, function(i){
    cseg.crop <- cseg==i
    r.crop <- EBImage::resize(r*cseg.crop, w=halfSize, h=halfSize)
    g.crop <- EBImage::resize(g*cseg.crop, w=halfSize, h=halfSize)
    b.crop <- EBImage::resize(b*cseg.crop, w=halfSize, h=halfSize)
    y.crop <- EBImage::resize(y*cseg.crop, w=halfSize, h=halfSize)
    EBImage::writeImage(r.crop, paste0(outHeader1, "_red_crop_", i, ".png"))
    EBImage::writeImage(g.crop, paste0(outHeader1, "_green_crop_", i, ".png"))
    EBImage::writeImage(b.crop, paste0(outHeader1, "_blue_crop_", i, ".png"))
    EBImage::writeImage(y.crop, paste0(outHeader1, "_yellow_crop_", i, ".png"))
  })
  
  ## A segmented RGB image for reference
  outHeader2 <- file.path(dataDir, "Segmented", cellID)
  seg <- EBImage::paintObjects(cseg, EBImage::paintObjects(nseg, img, col=c('#ffff00')), col=c('#ff00ff'))
  EBImage::writeImage(seg, paste0(outHeader2, "_segmented.png"))
  return(NULL)
}


