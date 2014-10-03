package pse;

public class VectorOfInt
{
    public VectorOfInt()
    {
        this(16);
    }

    public VectorOfInt(int capacity)
    {
        this.capacity = capacity;
        this.size = 0;
        this.data = new int[capacity];
    }

    public void pushBack(int val)
    {
        if (capacity <= size)
        {
            grow();
        }
        data[size++] = val;
    }

    public int capacity()
    {
        return capacity;
    }

    public int size()
    {
        return size;
    }

    public int[] data()
    {
        return data;
    }

    private void grow()
    {
        capacity *= 1.5;
        int[] dataOld = data;
        data = new int[capacity];
        System.arraycopy(dataOld, 0, data, 0, size);
    }

    private int capacity;
    private int size;
    private int[] data;
}
