import pya


def klayout_init_layout():
	# getting main references of the application
	try:
		app = pya.Application.instance()
	except Exception as e:
		return pya.Layout(), None

	mw = app.main_window()
	lv = mw.current_view()
	cv = None

	# this insures that lv and cv are valid objects
	if (lv == None):
		cv = mw.create_layout(1)
		lv = mw.current_view()
	else:
		cv = lv.active_cellview()

	# find or create the desired by programmer cell and layer
	return cv.layout(), lv


