/**********************************************************************
  Cp2kOutputDIalog - Analyze CP2k Output

  Copyright (C) Aoyama Iwao

  This file is not yet a part of the Avogadro molecular editor project.
  For more information about Avogadro, see <http://avogadro.cc/>

  Some code is based on Open Babel
  For more information, see <http://openbabel.sourceforge.net/>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 ***********************************************************************/

#include "cp2koutputdialog.h"

#include <QApplication>

namespace Avogadro
{
	#ifndef BOHR_TO_ANGSTROM
        #define BOHR_TO_ANGSTROM 0.529177249
    #endif


	Cp2kOutputDialog::Cp2kOutputDialog( QWidget * parent, Qt::WindowFlags flags ) 
		: QDialog( parent, flags ), m_molecule(0), m_animation(0), m_framesPerStep(8), m_vibScale(0.7), iMode(0)
	{

       ui.setupUi(this);

	   connect( ui.closeButton, SIGNAL(clicked()), this, SLOT(closeClicked()) );
	   connect( ui.loadFileButton, SIGNAL(clicked()), this, SLOT(loadFile()) );

	}

    Cp2kOutputDialog::~Cp2kOutputDialog()
	{
	}

	void Cp2kOutputDialog::setMolecule(Molecule *molecule)
    {
        m_molecule = molecule;
	}


	void Cp2kOutputDialog::createVibrationAnimation( )
	{
			 if (!m_animation)
                   m_animation = new Animation;
			 
			 else 
			 {
                   delete m_animation;
                   m_animation = new Animation;
				   
				   if(animationButton && pauseButton )
				   {
			           animationButton->disconnect();
			           pauseButton->disconnect();
				   }
             }

			 if( vibrationTable )
			 {
				 if( !(vibrationTable->selectedItems().isEmpty()) )
			       iMode = vibrationTable->selectedItems().at(0)->row();
			 }

			 unsigned int  nAtoms = m_molecule->numAtoms();
             QList<Atom *> atoms = m_molecule->atoms();

			 m_molecule->setConformer(0);  

			 // std::vector<std::vector<Eigen::Vector3d> *> m_curFrames;
			 m_curFrames.clear();

			 for(int i = 0; i< m_framesPerStep ; i++ )
				 m_curFrames.push_back(new std::vector<Eigen::Vector3d>(nAtoms));

             for(int i = 0; i< m_framesPerStep ; i++ )
			 {
				 for(int j=0; j<nAtoms; j++)
				 {
					 m_curFrames.at(i)->at(j).x() = atoms.at(j)->pos()->x() + vLx.at(iMode).at(j).x() * m_vibScale * sin(2*M_PI*i/(m_framesPerStep-1));
					 m_curFrames.at(i)->at(j).y() = atoms.at(j)->pos()->y() + vLx.at(iMode).at(j).y() * m_vibScale * sin(2*M_PI*i/(m_framesPerStep-1));
					 m_curFrames.at(i)->at(j).z() = atoms.at(j)->pos()->z() + vLx.at(iMode).at(j).z() * m_vibScale * sin(2*M_PI*i/(m_framesPerStep-1));
				 }
			 }

			  m_animation->setFrame(1); //Set the current frame
			  m_animation->setFps(10);
			  m_animation->setLoopCount(0); // argument 0 means "repeat forever".
			  m_animation->setFrames(m_curFrames);
			  m_animation->setMolecule(m_molecule);
	
			  if(animationButton && pauseButton )
		     {
			    connect(animationButton, SIGNAL(clicked()), m_animation, SLOT(start()) );
			    connect(pauseButton, SIGNAL(clicked()), m_animation, SLOT(pause()) );
			  }

			  // qDebug() << iMode << "th row is selected.";

			  return;

	}

	void Cp2kOutputDialog::closeClicked()
	{
		close();
	}

	void Cp2kOutputDialog::loadFile()
	{
			 // molecular geometry info
			 unsigned int nAtom = 0;
			 std::vector<std::string> a; // atoms
			 std::vector<vector3> coord; // coordinates

             // normal modes
			 unsigned int nNormMode = 0; // number of normal modes
			 std::vector< vector3 > curModeNormCoord; // atomic displacements of a normal mode

			  // Try to set "defaultpath" for dialog using the next sequence:
              //  (1) directory of current file (if any);
              //  (2) directory where previous deck was saved;
              //  (3) $HOME
		 
			 QFileInfo defaultFile( m_molecule->fileName() ); // directory of current file
             QString defaultPath = defaultFile.canonicalPath();

             QStringList filters;
             QString selectedFilter;

            if(m_savePath == "") // directory where previous deck was saved is empty;
		    {
                if (defaultPath.isEmpty())
			      { defaultPath = QDir::homePath(); } // $HOME
		    }

            else 
		    {  defaultPath = m_savePath; } // directory where previous deck was saved;

            if (m_saveFilter != "") 
			 { selectedFilter = m_saveFilter; }

			// For vibration output
            filters << tr("Molden vibration file") + " (*.mol)";

			QString defaultFileName = defaultPath + '/' + defaultFile.fileName();

            QString fileName = QFileDialog::getOpenFileName( this, tr("Read vibration file"),
                                                             defaultFileName, 
													         filters.join(";;"),
													         &selectedFilter);
			 if ( fileName == "" ) 
			 {
				  //QMessageBox msgBox;
                  //msgBox.setText(tr("No file found."));
                  //msgBox.exec();
				 qDebug() << tr("No output file found.") << endl;
                 return ;
             }

			 m_saveFilter = selectedFilter;

             QFile file(fileName);
             QTextStream in(&file);

             m_savePath = QFileInfo(file).absolutePath();

			// For molden format vibration output
            if (QFileInfo(file).completeSuffix().contains("mol")) 
		    {
               if(!file.open(QIODevice::ReadOnly | QIODevice::Text)) 
			   {
				   qDebug() << "Molden vibration file open failure." << endl;
				   return;
			   }
			  
			   bool freqOK = false;
               bool coordOK = false;
			   bool normCoordOK = false;
			   bool intOK = false;

			   while (!in.atEnd()) 
			   {
                    QString outputText = in.readLine(); // Reads one line of text from the stream

					// read frequencies
					if( outputText.contains("[FREQ]") )
					 { 
						 freqOK = true;
						 vF.clear();
						 outputText = in.readLine();
					 }

					double curFreq = 0.0;

					if( freqOK )
					{  curFreq = outputText.toDouble( &freqOK ); }

					if( freqOK )
					{ vF.push_back( curFreq ); }
					// end read frequencies

                    // read coordinate
					if( outputText.contains("[FR-COORD]") )
					{
						coordOK = true;
						outputText = in.readLine();
					}

					vector3 curCoord;
					QStringList infoText;
	
				    infoText = outputText.split(" ", QString::SkipEmptyParts);

					if( infoText.empty() || infoText.size() < 4 )
					 { coordOK = false;}

					if( coordOK )
					{
					  curCoord.SetX(infoText.at(1).toDouble( &coordOK ) * BOHR_TO_ANGSTROM ) ; 
					  curCoord.SetY(infoText.at(2).toDouble( &coordOK ) * BOHR_TO_ANGSTROM ) ; 
					  curCoord.SetZ(infoText.at(3).toDouble( &coordOK ) * BOHR_TO_ANGSTROM ) ; 
					}

					if( coordOK )
					{
						nAtom++;
						a.push_back(infoText.at(0).toStdString());
						coord.push_back(curCoord);
					}				
					// end read coordinate

					// read normal mode coordinates (atom displacements)
					vector3 curNormCoord;
					
					if( outputText.contains("[FR-NORM-COORD]") )
					{
						normCoordOK = true;
						vLx.clear();

						outputText = in.readLine(); // Skip "vibration 1" line.
						outputText = in.readLine();
					}

					if( outputText.contains("vibration") || outputText.contains("[INT]"))
					{ 
						vLx.push_back( curModeNormCoord );
						curModeNormCoord.clear();
						nNormMode++;
					}

					infoText = outputText.split(" ", QString::SkipEmptyParts);

					if( infoText.empty() || infoText.size() < 3 )
					{  normCoordOK = false; }

					if( normCoordOK )
					{
						curNormCoord.SetX( infoText.at(0).toDouble( &normCoordOK ) * BOHR_TO_ANGSTROM );
						curNormCoord.SetY( infoText.at(1).toDouble( &normCoordOK ) * BOHR_TO_ANGSTROM );
						curNormCoord.SetZ( infoText.at(2).toDouble( &normCoordOK ) * BOHR_TO_ANGSTROM );
					}

					if( normCoordOK )
					{	curModeNormCoord.push_back( curNormCoord );}

					if( outputText.contains("vibration") ) // For a next cycle.
					{ normCoordOK = true; }

					else if( outputText.contains("[INT]") ) // end [FR-NORM-COORD] block.
					{ normCoordOK = false;}
					// end read normal mode coordinates (atom displacements)

					// read intensities
					if( outputText.contains("[INT]") )
					{
						intOK = true;
						vI.clear();
						outputText = in.readLine();
					}

					double curInt = 0.0;

					if( intOK )
					{ curInt = outputText.toDouble( &intOK ); }

					if( intOK )
					{ vI.push_back( curInt ); }
					// end read intensities
			   }

			   // Set geometry found in vibration file(.mol) to m_molecule
			   OBMol* mol = new OBMol;
			
			   mol->BeginModify();

			   for( int i = 0; i < nAtom ; i++ )
			  {
				  OBAtom atom;
				  atom.SetAtomicNum( etab.GetAtomicNum( a.at(i).c_str() ));
				  atom.SetVector( coord.at(i) );
				  mol->AddAtom( atom );
			   }

			   mol->ConnectTheDots();
		       mol->PerceiveBondOrders(); 

			   mol->EndModify();

			   m_molecule->setOBMol(mol);
			   // end set geometry to m_molecule

			   if( !getVibrationWidgets() ) 
			   {
				   qDebug() << "Failed to get vibration widgets";
				   return;
			   }

			   // Show the vibration dock if it is hidden
			   if (!vibrationDock->isVisible()) 
			        vibrationDock->show();

			   // Activate the vibration widget
			   vibrationWidget->setEnabled(true);

			   // Activate animation buttons
			   animationButton->setEnabled(true);
			   pauseButton->setEnabled(true);

			   // Disconnect SIGNAL/SLOTs to avoid exception when a component is clicked.
			   vibrationTable->disconnect();
			   spectraButton->disconnect();
			   displayForcesCheckBox->disconnect();
			   animationButton->disconnect();
			   pauseButton->disconnect();

			   // Fill vibration table with vF and vI data.
			   vibrationTable->horizontalHeader()->show();
			   vibrationTable->setColumnCount(2);
			   vibrationTable->setRowCount(vF.size());

			   for (unsigned int row = 0; row < vF.size(); ++row) 
			   {
                    char buf[64];

                    sprintf(buf, "%lf", vF[row]);
                    QString Freq(buf);

					if( !vI.empty() )
                       sprintf(buf, "%lf", vI[row]);

					else sprintf(buf, "%lf", "0.0");

                    QString Int(buf);
	
                    QTableWidgetItem* newFreq = new QTableWidgetItem;
                    newFreq->setText(Freq);
                    QTableWidgetItem* newInten = new QTableWidgetItem;
                    newInten->setText(Int);

                    vibrationTable->setItem(row, 0, newFreq);
                    vibrationTable->setItem(row, 1, newInten);
			   }

			   // Temporary the following widgets are set disabled
			   editFilter->setEnabled(false);
			   filterLabel->setEnabled(false);
			   kmmolLabel->setEnabled(false);
			   scaleSlider->setEnabled(false);
		       normalizeDispCheckBox->setEnabled(false);
		       displayForcesCheckBox->setEnabled(false);
		       animationSpeedCheckBox->setEnabled(false);

			   // After table items are sorted by clicking a header, selection becomes invalid.
			   vibrationTable->horizontalHeader()->setEnabled(false);

			   createVibrationAnimation(); // For initial no-selection case

			   connect( vibrationTable, SIGNAL(itemSelectionChanged () ),
			   	        this, SLOT(createVibrationAnimation()) );

			   // Plot IR spectra
			   if( !vI.empty() )
			   {
			      SpectraPlotDialog *m_dialog = new SpectraPlotDialog(qobject_cast<QWidget*>(parent()));

			     if (m_molecule)
                     m_dialog->setIRSpectra(m_molecule, vF, vI);

			     connect( spectraButton, SIGNAL(clicked()), m_dialog, SLOT(show()) );
			   }


		   } // end Molden format vibration output read

		  // For debug
          /*
		   qDebug() << "vF:";
          for( int i = 0; i <vF.size(); i++ )
          {
			  // qDebug() << vF[i] << endl;
			  qDebug() << vF[i];
		  }
		  
		  qDebug() << "Atomic coordinates:"	 ;     
		  for( int i = 0; i < nAtom; i++ )
          {
			  qDebug() << a.at(i).c_str() << "  " << coord.at(i).x() << "  " << coord.at(i).y() << "  " << coord.at(i).z();
		  }
		  
		  qDebug() << "vLx:";
		   for(int i = 0; i<nNormMode; i++ )
		   {
			   qDebug() << i+1 << "th normal mode coordinates:";

			   for( int j = 0; j<nAtom ; j++ )
			   {
				   qDebug() << vLx.at(i).at(j).x() << "   " << vLx.at(i).at(j).y() << "   " << vLx.at(i).at(j).z();
			   }
		   }

		   qDebug() << "vI:";
	      for( int i = 0; i <vI.size(); i++ )
          {

			  qDebug() << vI[i];
		  }
          */
		  // end debug

			closeClicked();

			return;
			
	}// end loadFile()


	bool Cp2kOutputDialog::getVibrationWidgets()
	{
			   foreach (QWidget *widget, QApplication::allWidgets())
               {
                   if (widget->objectName().contains("vibrationDock",Qt::CaseInsensitive))
				   {
                        vibrationDock = widget;

                        vibrationWidget = vibrationDock->findChild<QWidget*>("VibrationWidget");
						  if(!vibrationWidget) return false;
                        vibrationTable = vibrationDock->findChild<QTableWidget*>("vibrationTable");
						  if(!vibrationTable) return false;
					    editFilter = vibrationDock->findChild<QLineEdit*>("editFilter");
						  if(!editFilter) return false;
						filterLabel = vibrationDock->findChild<QLabel*>("label_2");
						  if(!filterLabel) return false;
						kmmolLabel = vibrationDock->findChild<QLabel*>("label_3");
						  if(!kmmolLabel) return false;
                        spectraButton = vibrationDock->findChild<QPushButton*>("spectraButton");
						  if(!spectraButton) return false;
						scaleSlider = vibrationDock->findChild<QSlider*>("scaleSlider");
						  if(!scaleSlider) return false;
						normalizeDispCheckBox = vibrationDock->findChild<QCheckBox*>("normalizeDispCheckBox");
						  if(!normalizeDispCheckBox) return false;
                        displayForcesCheckBox = vibrationDock->findChild<QCheckBox*>("displayForcesCheckBox");
						  if(!displayForcesCheckBox) return false;
						animationSpeedCheckBox = vibrationDock->findChild<QCheckBox*>("animationSpeedCheckBox");
						  if(!animationSpeedCheckBox) return false;
                        animationButton = vibrationDock->findChild<QPushButton*>("animationButton");
						  if(!animationButton) return false;
                        pauseButton = vibrationDock->findChild<QPushButton*>("pauseButton");
						  if(!pauseButton) return false;
				   }

                }

			   if(!vibrationDock) return false;

			   else return true;

	}


   /*
   //----------------------------------------
   // SpectraPlotDialog class function definition
   //----------------------------------------
   */

   // Constructor
    SpectraPlotDialog::SpectraPlotDialog(QWidget *parent, Qt::WindowFlags f )
	    : QDialog( parent, f )
	{
		ui.setupUi(this);
	}

	SpectraPlotDialog::~SpectraPlotDialog() {};

	void SpectraPlotDialog::setIRSpectra(Molecule *mol, std::vector< double > vF, std::vector< double > vI)
	{
		ui.plot->resetPlot();

		m_molecule = mol;

        if (!m_molecule) return;
		if ( vI.empty() ) return;

		PlotObject *data = new PlotObject(Qt::red,  // the color for the plotting object
			                              PlotObject::Lines, // Plotting object. Points, Lines or Bars
			                                2); // the size to use for plotted points, in pixels 

	    double curX = vF.at(0);
		double minX = vF.at(0);
		double maxX = vF.at(0);

		double curY = vI.at(0);
		double minY = vI.at(0);
		double maxY = vI.at(0);

		// Find plotting boundaries
		for(int i = 0; i< vF.size() ; i++)
		{
			curX = vF.at(i);
			curY = vI.at(i);

			if (curX < minX) minX = curX;
            if (curX > maxX) maxX = curX;
			if (curY < minY) minY = curY;
            if (curY > maxY) maxY = curY;
		}

		curX = minX;
		double linewidth = 100.0;

		// calculate spectrum with specified bandwidth
		while (1)
		{
			curY = 0.0;

			for (int i = 0; i< vF.size(); i++)
			{
				curY += vI.at(i) * exp( -(1.0/linewidth) * (curX - vF.at(i)) * (curX - vF.at(i)));
			}

			data->addPoint(curX, -(curY), 0, 0.0); // set curY negative to turn the spectrum downward

			curX += 1.0; // increment a plot step
			if (curX > maxX) break;

		}


        double spreadX = maxX - minX;
		double spreadY = maxY - minY;
        double extX = spreadX * 0.05;
		double extY = spreadY * 0.05;

	   ui.plot->setDefaultLimits(minX - extX,
		                        maxX + extX,
		                        -( maxY + extY),
		                         -( minY - extY) );

	    ui.plot->setAntialiasing(true);
        ui.plot->setMouseTracking(true);
        ui.plot->axis(PlotWidget::BottomAxis)->setLabel(tr("Wavenumber(cm-1)"));
        ui.plot->axis(PlotWidget::LeftAxis)->setLabel(tr("-Intensity (km/mol)"));
		ui.plot->addPlotObject(data);
	}


}// end namespace Avogadro




