#include "Cell.h"

Cell::Cell()
{
    
}

Cell::~Cell()
{
    
}

void Cell::initCell(int i_index)
{
    m_index = i_index;
    m_headOfList = -9;
    m_numberOfNeighbors = 0;
}

void Cell::setHeadOfList(int i_headOfList)
{
    m_headOfList = i_headOfList;
}

void Cell::addNeighbor(int i_neighborIndex)
{
    m_listOfNeighbors[m_numberOfNeighbors] = i_neighborIndex;
    m_numberOfNeighbors++;
}

int Cell::headOfList()
{
    return m_headOfList;
}

int Cell::neighbor(int i_n)
{
    return m_listOfNeighbors[i_n];
}

int Cell::numberOfNeighbors()
{
    return m_numberOfNeighbors;
}
