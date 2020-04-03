class Cell
{
    
private:
    int m_index;
    int m_headOfList;
    int m_listOfNeighbors[9];
    int m_numberOfNeighbors;
    
public:
    Cell();
    ~Cell();
    void initCell(int);
    void setHeadOfList(int);
    void addNeighbor(int);
    
    int headOfList();
    int neighbor(int);
    int numberOfNeighbors();
};
