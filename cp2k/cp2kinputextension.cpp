/**********************************************************************
  Cp2kInputExtension - CP2k Input Extension

  Copyright (C) Aoyama Iwao

  This file is part of the Avogadro molecular editor project.
  For more information, see <http://avogadro.cc/>

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
#include "psfpotdialog.h"
#include "cp2kinputextension.h"

#include <QAction>
#include <QDebug>

using namespace std;

namespace Avogadro
{

  // this is a trick to identify what action we are taking
  enum Cp2kInputExtensionIndex
  {
    psfpot = 0,
    cp2kinput
  };

  Cp2kInputExtension::Cp2kInputExtension( QObject *parent ) : Extension( parent ) 
  {
    QAction* action;

	/*
	action = new QAction( this );
	action->setSeparator( true );
	m_actions.append( action );
	*/

	action = new QAction( this );
    action->setText( tr("&Psf/Pot..." ));
	action->setData( psfpot );
    m_actions.append( action );

	action = new QAction( this );
    action->setText( tr("&CP2K Input..." ));
	action->setData( cp2kinput );
    m_actions.append( action );

    action = new QAction( this );
	action->setSeparator( true );
	m_actions.append( action );

	// label = new QLabel("Would crate psf/pot!"); // for debug

  }

  Cp2kInputExtension::~Cp2kInputExtension()
  {
  }

  QList<QAction *> Cp2kInputExtension::actions() const
  {
    return m_actions;
  }

  // allows us to set the intended menu path for each action
  /*
  QString Cp2kInputExtension::menuPath(QAction *action) const
  {

    int i = action->data().toInt();

    switch ( i )
	{
      case psfpot:
        return tr("E&xtensions") + '>' + tr("&psf/pot...");
        break;
      case cp2kinput:
        return tr("E&xtensions") + '>' + tr("&cp2k input...");
        break;
    }
    return "";

  }
  */

  QDockWidget * Cp2kInputExtension::dockWidget()
  {
    // if we need a dock widget we can set one here
    return 0;
  }

  void Cp2kInputExtension::setMolecule(Molecule *molecule)
  {
    m_molecule = molecule;
  }

  QUndoCommand* Cp2kInputExtension::performAction(QAction *action, GLWidget *)
  {

    int i = action->data().toInt();

      PsfPotDialog* ppDialog = new PsfPotDialog(static_cast<QWidget*>(parent()));
 	  Cp2kInputDialog* cDialog = new Cp2kInputDialog(static_cast<QWidget*>(parent()));

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


		if (!cDialog) 
		{
           qDebug() << "No cp2k dialog ! Something went wrong!";
           return 0;
         }

		if(m_molecule)
             cDialog->setMolecule(m_molecule);

		cDialog->show();

        break;
    }


    return 0;
  }

}

Q_EXPORT_PLUGIN2(cp2kinputextension, Avogadro::Cp2kInputExtensionFactory)

