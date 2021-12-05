function test1(name, n)
    while n > 0
        println("$name $n")
        n -= 1
        yield()
    end
    false
end

t1 = @task test1("foo", 5)
t2 = @task test1("bar", 6)
t3 = @task test1("baz", 4)
schedule(t1)
schedule(t2)
schedule(t3)
wait(t2)