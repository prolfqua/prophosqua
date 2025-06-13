library(qrcode)

# 1) Generate the QR code object
doi <- "https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00441"
code <- qr_code(doi) # main qr_code() function :contentReference[oaicite:0]{index=0}

# 2) Plot with a caption
png("doi_qrcode_with_caption_prolfqua.png", width = 400, height = 450)
par(mar = c(5, 1, 1, 1)) # leave room at bottom
plot(code) # draw the QR code
mtext("prolfqua : 10.1021/acs.jproteome.2c00441",
  side = 1, line = 2, cex = 1.5
)
dev.off()


# 1) Generate the QR code object
doi <- "https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00911"
code <- qr_code(doi) # main qr_code() function :contentReference[oaicite:0]{index=0}

# 2) Plot with a caption
png("doi_qrcode_with_caption_prolfquapp.png", width = 400, height = 450)
par(mar = c(5, 1, 1, 1)) # leave room at bottom
plot(code) # draw the QR code
mtext("prolfquapp : 10.1021/acs.jproteome.4c00911",
  side = 1, line = 2, cex = 1.5
)
dev.off()
