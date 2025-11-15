progress_bar <- function(current,total)
{
	persentage = (current*100) / total
	current_width = (persentage*40)/100
	total_width = 40 - current_width
	bar <- paste0(strrep("=", current_width), strrep("-", total_width - current_width))
	sprintf("\r[%s] %3d%%", bar, persentage)
	flush.console()
}
