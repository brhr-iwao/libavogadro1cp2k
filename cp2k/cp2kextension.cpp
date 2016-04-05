/**********************************************************************
  Cp2kExtension - CP2k Extension

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

#include "cp2kinputdialog.h"
#include "cp2koutputdialog.h"
#include "psfpotdialog.h"
#include "Cp2kExtension.h"

#include <QAction>
#include <QDebug>
#include <QMessageBox> // To use QMessageBox for debug
#include <QDockWidget>

// Forward declaration of Avogadro::Molecule
/*
namespace Avogadro {
  class Molecule;
  class VibrationWidget;
}
*/

using namespace std;

namespace Avogadro
{

  // this is a trick to identify what action we are taking
  enum Cp2kExtensionIndex
  {
    psfpot = 0,
    cp2kinput,
	cp2koutput
  };

  Cp2kExtension::Cp2kExtension( QObject *parent ) : Extension( parent ) 
  {
    QAction* action;

	action = new QAction( this );
	action->setSeparator( true );
	m_actions.append( action );
	

	action = new QAction( this );
    action->setText( tr("&Psf/Pot..." ));
	action->setData( psfpot );
    m_actions.append( action );

	action = new QAction( this );
    action->setText( tr("Generate &Input..." ));
	action->setData( cp2kinput );
    m_actions.append( action );

    action = new QAction( this );
    action->setText( tr("Analyze &Output..." ));
	action->setData( cp2koutput );
    m_actions.append( action );

    action = new QAction( this );
	action->setSeparator( true );
	m_actions.append( action );

  }

  Cp2kExtension::~Cp2kExtension()
  {
  }

  QList<QAction *> Cp2kExtension::actions() const
  {
    return m_actions;
  }

  // allows us to set the intended menu path for each action
  QString Cp2kExtension::menuPath(QAction *action) const
  {

    int i = action->data().toInt();

    switch ( i )
	{
		/*
      case psfpot:
        return tr("E&xtensions") +'>'+ tr("&psf/pot...");
        break;
      case cp2kinput:
        return tr("E&xtensions") +'>'+tr("&CP2K") + '>' + tr("cp2k &input...");
        break;
	  case cp2koutput:
		return tr("E&xtensions") +'>'+tr("&CP2K") + '>' + tr("cp2k &output...");
        break;
		*/

      case psfpot:
        return tr("E&xtensions") ;
        break;
      case cp2kinput:
	  case cp2koutput:
        return tr("E&xtensions") +'>'+tr("&CP2K") ;
        break;

    }

    return "";

  }


  QDockWidget * Cp2kExtension::dockWidget()
  {
    // if we need a dock widget we can set one here
    return 0;
  }

  void Cp2kExtension::setMolecule(Molecule *molecule)
  {
    m_molecule = molecule;
  }

  QUndoCommand* Cp2kExtension::performAction(QAction *action, GLWidget *)
  {

    int i = action->data().toInt();

      PsfPotDialog* ppDialog = new PsfPotDialog(static_cast<QWidget*>(parent()));
 	  Cp2kInputDialog* ciDialog = new Cp2kInputDialog(static_cast<QWidget*>(parent()));
	  Cp2kOutputDialog* coDialog = new Cp2kOutputDialog(static_cast<QWidget*>(parent()));

    switch ( i )
	{
      case  psfpot:
        // perform action

		if (!ppDialog) 
		{
           qDebug() << "No Psf/Pot dialog ! Something went wrong!";
           return 0;
         }

		if(m_molecule)
             ppDialog->setMolecule(m_molecule);

		ppDialog->show();

        break;

      case cp2kinput:
        // perform action

		if (!ciDialog) 
		{
           qDebug() << "No cp2k input dialog ! Something went wrong!";
           return 0;
         }

		if(m_molecule)
             ciDialog->setMolecule(m_molecule);

		ciDialog->show();

        break;

	  case cp2koutput:

		if (!coDialog) 
		{
           qDebug() << "No cp2k outout dialog ! Something went wrong!";
           return 0;
         }

		if(m_molecule)
             coDialog->setMolecule(m_molecule);

		coDialog->show();

		break;
    }


    return 0;
  }

}

Q_EXPORT_PLUGIN2(Cp2kExtension, Avogadro::Cp2kExtensionFactory)

