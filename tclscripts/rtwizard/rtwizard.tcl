#
# Danube.tcl
#
# This is the main script for the RtWizard. 
#
# The development version codename is "Danube"
# The final version will be "Defiant" ;-)

#
# Extend Autopath
#
if { [info exists env(BRLCAD_ROOT)] } {
    set rootDir $env(BRLCAD_ROOT)
} else {
    set rootDir [file join / usr brlcad ]
} 
set wizDir [ file join $rootDir tclscripts rtwizard ]
set libDir [ file join $wizDir lib ]

lappend auto_path $wizDir $libDir
puts $auto_path

#
# All RtWizard stuff is in the RtWizard namespace
#
namespace eval RtWizard {

    #
    # Required packages. There are several!
    #
    package require Itcl
    package require Itk
    package require Iwidgets 3.0

    package require Wizard 1.0
    package require DbPage
    package require FbPage
    package require HelpPage
    package require ExamplePage

    package require FullColorPage
    package require HighlightedPage
    package require GhostPage
    package require LinePage

    package require FeedbackDialog

    #
    # The supported picture types
    #
    package require PictureTypeA
    package require PictureTypeB
    package require PictureTypeC
    package require PictureTypeD
    package require PictureTypeE
    package require PictureTypeF

    set dbFile ""
    set helpFont {-family helvetica -size 12 \
	    -slant italic}

    #
    # main
    #
    # Do all that is doable
    #
    proc main { args } {

	#
	# Create the Feedback
	#
	set fb [ FeedbackDialog .#auto \
		-title "Loading RtWizard ..." \
		-width 250 \
		-height 100 ]

	$fb center
	$fb activate
	$fb center

	#
	# Create the wizard
	#
	set w [ RtWizard::Wizard .#auto ]
	$fb inform "Wizard megawidget created" 5

	#
	# Add the first set of pages to the wizard.
	# Note that we can't create the Mged-based pages until
	# the database is selected.
	#
	# XXX Fixing this will require debugging the Mged object
	#     The opendb command has a fatal bug.
	#

	# Database Page
	set ::dbp [$w add RtWizard::DbPage "dbp"]
	$w enable "dbp"
	$fb inform "Database page created." 6

	# Framebuffer Page
	set ::fbp [$w add RtWizard::FbPage "fbp"]
	$w enable "fbp"
	$fb inform "Frame buffer page created." 4

	$w add RtWizard::IntroPage "intro"
	$w enable "intro"
	$fb inform "Intro page created." 4

	$w add RtWizard::HelpPage "help"
	$w enable "help"
	$fb inform "Help page created" 4

	set ::exp [$w add RtWizard::ExamplePage "exp"]
	$w enable "exp"
	$fb inform "Example page created" 4

	#
	# Load the picture types into the example page.
	# It then forwards the info to the wizard for the 
	# radiobox in the "Image" menu. The "switch" for
	# the picture type is maintained by the ExamplePage
	# 
	$::exp addType RtWizard::PictureTypeA typeA
	$fb inform "Support for image type A loaded." 1

	$::exp addType RtWizard::PictureTypeB typeB
	$fb inform "Support for image type B loaded." 1

	$::exp addType RtWizard::PictureTypeC typeC
	$fb inform "Support for image type C loaded." 1

	$::exp addType RtWizard::PictureTypeD typeD
	$fb inform "Support for image type D loaded." 1

	$::exp addType RtWizard::PictureTypeE typeE
	$fb inform "Support for image type E loaded." 2

	$::exp addType RtWizard::PictureTypeF typeF
	$fb inform "Support for image type F loaded." 2

	#
	# XXX need to disable all menu options 
	#

	#
	# Pack the wizard
	#
	pack $w
	$fb inform "Megawidgets packed." 5

	#
	# Get the database name. If the name was specified on the 
	# command line, use it and proceed. Otherwise, spin up the 
	# gui, select the database page, and wait for the database
	# file to be specified.
	#
	if { [llength $args] > 0 } {
	    set ::RtWizard::dbFile [ lindex $args 0 ]
	    if { ! [file exists $::RtWizard::dbFile] } {
		set ::RtWizard::dbFile ""
	    }
	} 
	
	if { [string length $::RtWizard::dbFile] == 0 } {
	    #
	    # select the database page
	    # 
	    $w select "dbp"

	    #
	    # Hide the Feedback
	    #
	    wm iconify $fb

	    #
	    # Start up the gui, and run until the dbFile has
	    # be specified.
	    #
	    vwait ::RtWizard::dbFile

	    #
	    # Restore the Feedback
	    #
	    wm deiconify $fb
	}

	#
	# At this point, the database should be specified. 
	# We load the database ourselves and hand it to the 
	# pages.
	#
	set ::mgedObj [Mged .#auto $::RtWizard::dbFile]
	$::mgedObj configure -multi_pane 0
	$fb inform "MGED object instantiated." 40

	$w add RtWizard::FullColorPage fullColor $::RtWizard::dbFile
	$fb inform "Support for full color images loaded." 5

	$w add RtWizard::HighlightedPage highlighted $::RtWizard::dbFile
	$fb inform "Support for highlighted images loaded." 5

	$w add RtWizard::GhostPage ghost $::RtWizard::dbFile
	$fb inform "Support for ghost images loaded." 5

	$w add RtWizard::LinePage lines $::RtWizard::dbFile 
	$w select "exp"
	$fb inform "RtWizard ready!!." 5

	#
	# Pause...
	#
	after 1000
	$fb deactivate

	#
	# Go!
	#
	vwait forever
    }

}

#
# Start main
# 
RtWizard::main $argv














